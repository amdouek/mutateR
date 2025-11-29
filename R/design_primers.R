#' Design Genotyping Primers for Exon Deletions
#'
#' Scans genomic regions flanking the predicted deletion site to design
#' Forward and Reverse PCR primers. Uses a penalty-based selection system
#' optimizing for Tm (~60°C), GC content (50%), and 3' GC clamps.
#'
#' @param pairs_df Data frame returned by \code{assemble_grna_pairs}.
#' @param genome BSgenome object.
#' @param flank_min Integer. Minimum distance from cut site to primer (default 100bp).
#'                  Ensures primers are not chewed back during repair.
#' @param flank_max Integer. Maximum distance from cut site to primer (default 400bp).
#' @param target_tm Numeric. Optimal melting temperature (default 60).
#'
#' @return The input data frame with 4 new columns:
#'   \item{primer_fwd}{Sequence of the forward primer (5'-3')}
#'   \item{primer_rev}{Sequence of the reverse primer (5'-3')}
#'   \item{wt_amplicon_size}{Size of the band in wild-type cells}
#'   \item{mut_amplicon_size}{Predicted size of the band in deletion clones}
#'
#' @importFrom Biostrings getSeq reverseComplement DNAString matchPattern
#' @export
get_genotyping_primers <- function(pairs_df,
                                   genome,
                                   flank_min = 100,
                                   flank_max = 400,
                                   target_tm = 60.0) {

  if (is.null(pairs_df) || nrow(pairs_df) == 0) return(pairs_df)
  if (!inherits(genome, "BSgenome")) stop("genome must be a BSgenome object.")

  message("Designing genotyping primers for ", nrow(pairs_df), " candidate pairs...")

  # --- Pre-allocate result vectors ---
  fwd_seqs <- character(nrow(pairs_df))
  rev_seqs <- character(nrow(pairs_df))
  wt_sizes <- integer(nrow(pairs_df))
  mut_sizes <- integer(nrow(pairs_df))

  # --- Helper: Simple Tm Calculation (GC-adjusted) ---
  # Formula: Tm = 64.9 + 41 * (GC - 16.4) / L
  calculate_tm <- function(seq_str) {
    len <- nchar(seq_str)
    if (len == 0) return(0)
    g_c <- sum(charToRaw(seq_str) %in% charToRaw("GC"))
    tm <- 64.9 + 41 * ((g_c - 16.4) / len)
    return(tm)
  }

  # --- Helper: Scan a region for the best primer ---
  # direction 1 = Forward (Genomic + strand)
  # direction -1 = Reverse (Genomic - strand, output as RevComp)
  find_best_primer <- function(chrom, start, end, direction = 1) {

    # 1. Get Sequence
    # We always fetch the + strand genomic sequence first
    region_seq <- tryCatch({
      as.character(Biostrings::getSeq(genome, names=chrom, start=max(1, start), end=end))
    }, error = function(e) return(NULL))

    if (is.null(region_seq) || nchar(region_seq) < 20) return(NULL)

    # 2. Generate Candidates (Sliding Window)
    # Lengths 18 to 24 bp
    candidates <- list()

    for (len in 18:24) {
      if (nchar(region_seq) < len) next
      num_windows <- nchar(region_seq) - len + 1
      starts <- 1:num_windows

      substrings <- substring(region_seq, starts, starts + len - 1)

      # 3. Filter & Score
      for (k in seq_along(substrings)) {
        s <- substrings[k]

        # GC Check (40-60%)
        gc_count <- sum(charToRaw(s) %in% charToRaw("GC"))
        gc_perc <- (gc_count / len) * 100
        if (gc_perc < 40 || gc_perc > 65) next

        # GC Clamp (Ends in G or C)
        last_base <- substring(s, len, len)
        if (!last_base %in% c("G", "C")) next

        # Poly-X Check (Max 4 repeats)
        if (grepl("AAAAA|CCCCC|GGGGG|TTTTT", s)) next

        # Tm Check
        tm <- calculate_tm(s)
        if (tm < 55 || tm > 65) next

        # Calculate Penalty
        # 1. Deviation from target Tm
        # 2. Deviation from 50% GC
        # 3. Distance from 'start' (we prefer primers closer to cut site to keep amplicon small)
        #    If Forward: 'k' increases = further from cut? No, region is [Cut-Max, Cut-Min].
        #    So higher 'k' = closer to cut. We want higher k.
        #    If Reverse: region is [Cut+Min, Cut+Max]. Lower 'k' = closer to cut.

        pen_tm <- abs(tm - target_tm) * 2
        pen_gc <- abs(gc_perc - 50) * 0.5

        candidates[[length(candidates) + 1]] <- list(
          seq = s,
          pos_relative = k,
          len = len,
          penalty = pen_tm + pen_gc
        )
      }
    }

    if (length(candidates) == 0) return(NULL)

    # Convert to DF
    cand_df <- do.call(rbind, lapply(candidates, as.data.frame))

    # 4. Pick Best
    # We want primers closer to the cut sites (smaller amplicons usually better for genotyping)
    # Forward Region: [Far ... Near] -> We want high relative position (end of string)
    # Reverse Region: [Near ... Far] -> We want low relative position (start of string)

    if (direction == 1) {
      # Forward: Prioritize higher 'pos_relative' (closer to 3' end of the search window)
      cand_df$final_score <- cand_df$penalty - (cand_df$pos_relative * 0.05)
    } else {
      # Reverse: Prioritize lower 'pos_relative' (closer to 5' end of the search window)
      cand_df$final_score <- cand_df$penalty + (cand_df$pos_relative * 0.05)
    }

    best <- cand_df[which.min(cand_df$final_score), ]

    # 5. Return info
    # Return genomic coordinate of the 5' end of the primer
    genomic_start <- start + best$pos_relative - 1

    res <- list(seq = as.character(best$seq),
                gen_start = genomic_start,
                gen_end = genomic_start + best$len - 1)

    if (direction == -1) {
      # Reverse Complement for Reverse Primer
      res$seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(res$seq)))
    }

    return(res)
  }

  # --- Main Loop ---
  for (i in seq_len(nrow(pairs_df))) {

    # Determine Cut Sites
    # Use explicit columns if available, or derive from start/end
    # Note: assemble_grna_pairs usually outputs 'cut_site_5p' / 'cut_site_3p' in newer versions,
    # but might fallback to 'end_5p' / 'start_3p' if simplified.

    # Logic: Cas9 cut is usually -3 from PAM.
    # Let's rely on genomic coords.
    # 5p site: We cut roughly at the end of the 5p guide.
    # 3p site: We cut roughly at the start of the 3p guide.

    chr <- pairs_df$seqnames_5p[i] # Assumes both on same chromosome

    # Define "Cut" coordinates roughly for search windows
    # (Exact cut site isn't strictly necessary for primer placement, just need to be flanking)

    if ("cut_site_5p" %in% names(pairs_df)) {
      c5 <- pairs_df$cut_site_5p[i]
      c3 <- pairs_df$cut_site_3p[i]
    } else {
      # Fallback approximation
      c5 <- pairs_df$end_5p[i]
      c3 <- pairs_df$start_3p[i]
    }

    # Ensure ordered (c5 < c3)
    start_cut <- min(c5, c3)
    end_cut   <- max(c5, c3)

    # --- Design Forward ---
    # Look upstream of start_cut
    # Region: [start_cut - flank_max, start_cut - flank_min]
    fwd_res <- find_best_primer(chr,
                                start = start_cut - flank_max,
                                end = start_cut - flank_min,
                                direction = 1)

    # --- Design Reverse ---
    # Look downstream of end_cut
    # Region: [end_cut + flank_min, end_cut + flank_max]
    rev_res <- find_best_primer(chr,
                                start = end_cut + flank_min,
                                end = end_cut + flank_max,
                                direction = -1)

    if (!is.null(fwd_res) && !is.null(rev_res)) {
      fwd_seqs[i] <- fwd_res$seq
      rev_seqs[i] <- rev_res$seq

      # WT Amplicon: Distance from Fwd 5' to Rev 5' (on genomic scale) + RevLen
      # Physically: (Rev_Gen_End) - (Fwd_Gen_Start) + 1
      # Note: Rev_res returns genomic coordinates of the binding site on the + strand.
      # Fwd is upstream (lower coord), Rev is downstream (higher coord).

      wt_size <- (rev_res$gen_end - fwd_res$gen_start) + 1
      wt_sizes[i] <- wt_size

      # Mutant Amplicon: WT - Deletion Size
      del_size <- if("genomic_deletion_size" %in% names(pairs_df)) pairs_df$genomic_deletion_size[i] else (end_cut - start_cut)
      mut_sizes[i] <- wt_size - del_size
    } else {
      fwd_seqs[i] <- NA_character_
      rev_seqs[i] <- NA_character_
      wt_sizes[i] <- NA_integer_
      mut_sizes[i] <- NA_integer_
    }
  }

  # --- Bind results ---
  pairs_df$primer_fwd <- fwd_seqs
  pairs_df$primer_rev <- rev_seqs
  pairs_df$wt_amplicon_size <- wt_sizes
  pairs_df$mut_amplicon_size <- mut_sizes

  n_success <- sum(!is.na(pairs_df$primer_fwd))
  message("Primers successfully designed for ", n_success, " / ", nrow(pairs_df), " pairs.")

  return(pairs_df)
}
