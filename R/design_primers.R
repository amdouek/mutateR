#' Design genotyping primers for mutateR-designed deletions
#'
#' Automatically designs PCR primers to genotype CRISPR-mediated deletions.
#' Uses the primer3-py backend for thermodynamic accuracy and strict filtering
#' (homopolymers, GC clamps, self-dimers).
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

  # --- Cache for Internal Primers (Strategy B) ---
  internal_primer_cache <- list()
  exons_df <- as.data.frame(exons_gr)

  message("Designing primers via Primer3 (Python backend)...")
  pb <- utils::txtProgressBar(min = 0, max = n_total, style = 3)

  for (i in seq_len(n_total)) {
    utils::setTxtProgressBar(pb, i)
    if (!has_coords[i]) next

    chrom <- as.character(pairs_df$seqnames_5p[i])
    c5 <- pairs_df$cut_site_5p[i]; c3 <- pairs_df$cut_site_3p[i]
    cut_s <- min(c5, c3); cut_e <- max(c5, c3)
    del_size <- abs(c3 - c5)

    # Padding for template extraction
    padding <- 600

    # Determine Strategy
    use_strat_b <- del_size > max_wt_amplicon

    # ---- Strategy A: Flanking Pair (Standard) ----
    if (!use_strat_b) {
      strategies[i] <- "Flanking"

      # 1. Fetch Template (WT)
      tmpl_start <- max(1, cut_s - padding)
      tmpl_end   <- cut_e + padding

      seq_dna <- tryCatch({
        as.character(Biostrings::getSeq(genome, names=chrom, start=tmpl_start, end=tmpl_end))
      }, error = function(e) NULL)

      if (!is.null(seq_dna) && nchar(seq_dna) > 0) {
        # 2. Define Target (Relative)
        # Target is the deletion span. Primer3 must design OUTSIDE this.
        rel_start <- cut_s - tmpl_start
        rel_len   <- cut_e - cut_s

        # 3. Call Python
        p_res <- run_primer3_python(
          sequence_template = seq_dna,
          target_start = rel_start,
          target_len = rel_len,
          tm_opt = target_tm,
          prod_min = rel_len + 80,   # Min amplicon = deletion + flanks
          prod_max = rel_len + 1000  # Max amplicon
        )

        if (!is.null(p_res)) {
          ext_fwd[i] <- p_res$fwd_seq
          ext_rev[i] <- p_res$rev_seq
          size_wt[i] <- as.character(p_res$prod_size)
          size_mut[i] <- p_res$prod_size - del_size
        }
      }

    } else {
      # ---- Strategy B: Dual Pair (Internal + External) ----
      strategies[i] <- "Dual_Pair"
      size_wt[i] <- "Too Large"

      # --- B1. External Pair (Mutant Band Detection) ---
      # Simulate the mutant allele by joining upstream and downstream sequences.
      # This allows Primer3 to "see" the post-deletion product size accurately.

      up_seq <- tryCatch(as.character(Biostrings::getSeq(genome, names=chrom, start=cut_s-400, end=cut_s)), error=function(e) "")
      dn_seq <- tryCatch(as.character(Biostrings::getSeq(genome, names=chrom, start=cut_e, end=cut_e+400)), error=function(e) "")

      if (nchar(up_seq) > 0 && nchar(dn_seq) > 0) {
        mut_templ <- paste0(up_seq, dn_seq)
        # Junction is at index 400 (length of up_seq)

        p_res <- run_primer3_python(
          sequence_template = mut_templ,
          target_start = 400, # Target the junction
          target_len = 0,     # Zero length target forces primers around this point
          tm_opt = target_tm,
          prod_min = 150,     # Standard small amplicon for genotyping
          prod_max = 600
        )

        if (!is.null(p_res)) {
          ext_fwd[i] <- p_res$fwd_seq
          ext_rev[i] <- p_res$rev_seq
          size_mut[i] <- p_res$prod_size
        }
      }

      # --- B2. Internal Pair (WT Band Detection) ---
      # Detects a region that SHOULD be deleted.

      e5_rank <- pairs_df$exon_5p[i]; e3_rank <- pairs_df$exon_3p[i]
      cache_key <- paste(e5_rank, e3_rank, sep="_")

      if (!is.null(internal_primer_cache[[cache_key]])) {
        cached <- internal_primer_cache[[cache_key]]
        int_fwd[i]  <- cached$fwd
        int_rev[i]  <- cached$rev
        size_int[i] <- cached$size
      } else {
        # Define internal region
        del_exons <- exons_df[exons_df$rank > e5_rank & exons_df$rank < e3_rank, ]

        target_s <- NA; target_e <- NA

        # Priority 1: Largest deleted exon (functional relevance)
        if (nrow(del_exons) > 0) {
          best_ex <- del_exons[which.max(del_exons$exon_chrom_end - del_exons$exon_chrom_start), ]
          w <- best_ex$exon_chrom_end - best_ex$exon_chrom_start
          if (w > 80) {
            target_s <- best_ex$exon_chrom_start
            target_e <- best_ex$exon_chrom_end
          }
        }

        # Priority 2: Midpoint of deletion (Intron)
        if (is.na(target_s)) {
          mid <- floor((cut_s + cut_e) / 2)
          target_s <- mid - 200
          target_e <- mid + 200
        }

        # Fetch sequence
        int_seq <- tryCatch({
          as.character(Biostrings::getSeq(genome, names=chrom, start=target_s, end=target_e))
        }, error = function(e) NULL)

        if (!is.null(int_seq) && nchar(int_seq) > 100) {
          # Target the middle of this sequence to ensure robust design
          mid_idx <- floor(nchar(int_seq) / 2)

          i_res <- run_primer3_python(
            sequence_template = int_seq,
            target_start = mid_idx,
            target_len = 0,
            tm_opt = target_tm,
            prod_min = 80,
            prod_max = min(400, nchar(int_seq))
          )

          # Cache and Assign
          if (!is.null(i_res)) {
            internal_primer_cache[[cache_key]] <- list(fwd=i_res$fwd_seq, rev=i_res$rev_seq, size=i_res$prod_size)
            int_fwd[i]  <- i_res$fwd_seq
            int_rev[i]  <- i_res$rev_seq
            size_int[i] <- i_res$prod_size
          } else {
            internal_primer_cache[[cache_key]] <- list(fwd=NA, rev=NA, size=NA)
          }
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
