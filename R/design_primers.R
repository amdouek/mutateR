#' Design genotyping primers for mutateR-designed deletions
#'
#' Automatically designs PCR primers to genotype CRISPR-mediated deletions.
#' Implements a dynamic strategy based on deletion size with caching for performance.
#'
#' @param pairs_df Data frame returned by \code{assemble_grna_pairs}.
#' @param exons_gr GRanges object of exon structures.
#' @param genome BSgenome object.
#' @param max_wt_amplicon Integer. Deletions larger than this trigger Strategy B (default 3000).
#' @param target_tm Numeric. Optimal melting temperature (default 60.0).
#'
#' @return The input data frame with added primer columns.
#' @export
get_genotyping_primers <- function(pairs_df,
                                   exons_gr,
                                   genome,
                                   max_wt_amplicon = 3000,
                                   target_tm = 60.0) {

  if (is.null(pairs_df) || nrow(pairs_df) == 0) return(pairs_df)

  # --- Input Validation ---
  if (!all(c("cut_site_5p", "cut_site_3p") %in% names(pairs_df))) {
    if ("end_5p" %in% names(pairs_df)) pairs_df$cut_site_5p <- pairs_df$end_5p
    if ("start_3p" %in% names(pairs_df)) pairs_df$cut_site_3p <- pairs_df$start_3p
  }

  has_coords <- !is.na(pairs_df$cut_site_5p) &
    !is.na(pairs_df$cut_site_3p) &
    !is.na(pairs_df$seqnames_5p)

  if (sum(!has_coords) > 0) warning("Skipping ", sum(!has_coords), " pairs due to missing coordinate data.")

  n_total <- nrow(pairs_df)

  # --- Setup Output Vectors ---
  strategies <- rep(NA_character_, n_total)
  ext_fwd    <- rep(NA_character_, n_total); ext_rev <- rep(NA_character_, n_total)
  int_fwd    <- rep(NA_character_, n_total); int_rev <- rep(NA_character_, n_total)
  size_wt    <- rep(NA_character_, n_total); size_mut <- rep(NA_integer_, n_total)
  size_int   <- rep(NA_integer_, n_total)

  # --- Cache for Internal Primers ---
  internal_primer_cache <- list()

  # --- Helper: Fast Tm Calc ---
  calc_tm_vec <- function(seqs) {
    w <- nchar(seqs)
    gc_counts <- nchar(gsub("[AT]", "", seqs))
    64.9 + 41 * ((gc_counts - 16.4) / w)
  }

  # --- Helper: Optimised Region Scanner ---
  scan_region <- function(chrom, start, end, direction = 1) {
    if (is.na(start) || is.na(end) || end <= start) return(NULL)

    seq_dna <- tryCatch({
      Biostrings::getSeq(genome, names=chrom, start=max(1, start), end=end)
    }, error = function(e) NULL)

    if (is.null(seq_dna) || length(seq_dna) == 0) return(NULL)

    seq_str <- as.character(seq_dna)
    len_region <- nchar(seq_str)
    p_len <- 20

    if (len_region < p_len) return(NULL)

    starts <- 1:(len_region - p_len + 1)
    cands  <- substring(seq_str, starts, starts + p_len - 1)

    # Filter GC (40-65%)
    gc_nums <- nchar(gsub("[AT]", "", cands))
    gc_frac <- gc_nums / p_len
    keep_gc <- gc_frac >= 0.40 & gc_frac <= 0.65
    if (!any(keep_gc)) return(NULL)
    cands_sub <- cands[keep_gc]; starts_sub <- starts[keep_gc]

    # Filter 3' Clamp
    last_char <- substring(cands_sub, p_len, p_len)
    keep_clamp <- last_char %in% c("G", "C")
    if (!any(keep_clamp)) return(NULL)
    cands_sub <- cands_sub[keep_clamp]; starts_sub <- starts_sub[keep_clamp]

    # Filter Tm
    tms <- calc_tm_vec(cands_sub)
    tm_diff <- abs(tms - target_tm)
    keep_tm <- tm_diff <= 5
    if (!any(keep_tm)) return(NULL)

    cands_final <- cands_sub[keep_tm]
    starts_final <- starts_sub[keep_tm]
    scores_final <- tm_diff[keep_tm]

    # Pick Best
    pos_norm <- starts_final / len_region
    pen_pos  <- if (direction == 1) (1 - pos_norm) * 2 else pos_norm * 2
    best_idx <- which.min(scores_final + pen_pos)

    final_seq <- cands_final[best_idx]
    final_gen_start <- start + starts_final[best_idx] - 1

    if (direction == -1) {
      final_seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(final_seq)))
    }
    return(list(seq = final_seq, start = final_gen_start, end = final_gen_start + p_len - 1))
  }

  exons_df <- as.data.frame(exons_gr)

  message("Designing primers (optimised)...")
  pb <- utils::txtProgressBar(min = 0, max = n_total, style = 3)

  for (i in seq_len(n_total)) {
    utils::setTxtProgressBar(pb, i)
    if (!has_coords[i]) next

    chrom <- as.character(pairs_df$seqnames_5p[i])
    c5 <- pairs_df$cut_site_5p[i]; c3 <- pairs_df$cut_site_3p[i]
    cut_s <- min(c5, c3); cut_e <- max(c5, c3)
    del_size <- abs(c3 - c5)

    # Determine Strategy
    use_strat_b <- del_size > max_wt_amplicon

    if (!use_strat_b) {
      # === Strategy A (Flanking) ===
      strategies[i] <- "Flanking"
      fwd <- scan_region(chrom, cut_s - 500, cut_s - 100, 1)
      rev <- scan_region(chrom, cut_e + 100, cut_e + 500, -1)
      if (!is.null(fwd) && !is.null(rev)) {
        ext_fwd[i] <- fwd$seq; ext_rev[i] <- rev$seq
        wt_len <- (rev$end - fwd$start + 1)
        size_wt[i]  <- as.character(wt_len)
        size_mut[i] <- wt_len - del_size
      }
    } else {
      # === Strategy B (Dual Pair) ===
      strategies[i] <- "Dual_Pair"
      size_wt[i] <- "Too Large"

      # 1. External Pair
      fwd <- scan_region(chrom, cut_s - 300, cut_s - 100, 1)
      rev <- scan_region(chrom, cut_e + 100, cut_e + 300, -1)
      if (!is.null(fwd) && !is.null(rev)) {
        ext_fwd[i] <- fwd$seq; ext_rev[i] <- rev$seq
        size_mut[i] <- (cut_s - fwd$start) + (rev$end - cut_e)
      }

      # 2. Internal Pair (Memoized)
      e5_rank <- pairs_df$exon_5p[i]; e3_rank <- pairs_df$exon_3p[i]
      cache_key <- paste(e5_rank, e3_rank, sep="_")

      if (!is.null(internal_primer_cache[[cache_key]])) {
        cached <- internal_primer_cache[[cache_key]]
        int_fwd[i] <- cached$fwd; int_rev[i] <- cached$rev; size_int[i] <- cached$size
      } else {
        # Determine internal target window
        del_exons <- exons_df[exons_df$rank > e5_rank & exons_df$rank < e3_rank, ]

        search_s <- NA; search_e <- NA
        use_midpoint <- TRUE

        # Try to find a valid deleted exon > 100bp
        if (nrow(del_exons) > 0) {
          best_ex <- del_exons[which.max(del_exons$end - del_exons$start), ]
          if ((best_ex$end - best_ex$start) > 100) {
            search_s <- best_ex$start
            search_e <- best_ex$end
            use_midpoint <- FALSE
          }
        }

        # Fallback to deletion midpoint (intron) if no valid exon
        if (use_midpoint) {
          mid <- floor((cut_s + cut_e) / 2)
          search_s <- mid - 200
          search_e <- mid + 200
        }

        # Perform Scan
        int_res <- NULL
        if (!is.na(search_s) && (search_e - search_s > 100)) {
          center <- floor((search_s + search_e) / 2)
          i_fwd <- scan_region(chrom, max(search_s, center-200), center-50, 1)
          i_rev <- scan_region(chrom, center+50, min(search_e, center+200), -1)
          if (!is.null(i_fwd) && !is.null(i_rev)) {
            int_res <- list(fwd = i_fwd$seq, rev = i_rev$seq, size = i_rev$end - i_fwd$start + 1)
          }
        }

        # Cache result
        internal_primer_cache[[cache_key]] <- if(is.null(int_res)) list(fwd=NA, rev=NA, size=NA) else int_res

        if (!is.null(int_res)) {
          int_fwd[i] <- int_res$fwd; int_rev[i] <- int_res$rev; size_int[i] <- int_res$size
        }
      }
    }
  }
  close(pb)

  # --- Assemble Output ---
  pairs_df$priming_strategy <- strategies
  pairs_df$primer_ext_fwd   <- ext_fwd
  pairs_df$primer_ext_rev   <- ext_rev
  pairs_df$primer_int_fwd   <- int_fwd
  pairs_df$primer_int_rev   <- int_rev
  pairs_df$exp_wt_size      <- size_wt
  pairs_df$exp_mut_size     <- size_mut
  pairs_df$exp_int_size     <- size_int

  n_succ <- sum(!is.na(pairs_df$primer_ext_fwd))
  message("\nDesigned primers for ", n_succ, "/", n_total, " pairs.")

  return(pairs_df)
}
