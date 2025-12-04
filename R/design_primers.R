#' Design genotyping primers for mutateR-designed deletions (Batched)
#'
#' Automatically designs PCR primers to genotype CRISPR-mediated deletions.
#' Uses the primer3-py backend for thermodynamic accuracy.
#' Optimized with batch processing to reduce R-Python overhead.
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

  # --- 1. Coordinate Validation ---
  if (!all(c("cut_site_5p", "cut_site_3p") %in% names(pairs_df))) {
    if ("end_5p" %in% names(pairs_df)) pairs_df$cut_site_5p <- pairs_df$end_5p
    if ("start_3p" %in% names(pairs_df)) pairs_df$cut_site_3p <- pairs_df$start_3p
  }

  has_coords <- !is.na(pairs_df$cut_site_5p) &
    !is.na(pairs_df$cut_site_3p) &
    !is.na(pairs_df$seqnames_5p)

  if (sum(!has_coords) > 0) {
    warning("Dropping ", sum(!has_coords), " pairs due to missing coordinate data.")
    pairs_df <- pairs_df[has_coords, ]
  }

  if (nrow(pairs_df) == 0) return(pairs_df)

  # --- 2. Filter Identical Cut Sites (0 bp Deletions) ---
  is_distinct <- pairs_df$cut_site_5p != pairs_df$cut_site_3p
  if (any(!is_distinct)) {
    pairs_df <- pairs_df[is_distinct, ]
  }

  if (nrow(pairs_df) == 0) return(pairs_df)

  # --- Setup Variables ---
  n_total <- nrow(pairs_df)
  exons_df <- as.data.frame(exons_gr)
  padding <- 600

  batch_requests <- list()
  internal_cache <- list()
  pending_assignments <- list()

  strategies <- rep(NA_character_, n_total)
  ext_fwd    <- rep(NA_character_, n_total); ext_rev <- rep(NA_character_, n_total)
  int_fwd    <- rep(NA_character_, n_total); int_rev <- rep(NA_character_, n_total)
  size_wt    <- rep(NA_character_, n_total); size_mut <- rep(NA_integer_, n_total)
  size_int   <- rep(NA_integer_, n_total)

  # --- Helper: Safe Sequence Fetching ---
  safe_get_seq <- function(...) {
    val <- tryCatch(as.character(Biostrings::getSeq(...)), error = function(e) "")
    if (length(val) == 0) return("")
    return(val)
  }

  message("Preparing primer design batch requests...")
  pb <- utils::txtProgressBar(min = 0, max = n_total, style = 3)

  # --- PHASE 1: Build Requests in R ---
  for (i in seq_len(n_total)) {
    utils::setTxtProgressBar(pb, i)

    chrom <- as.character(pairs_df$seqnames_5p[i])
    c5 <- pairs_df$cut_site_5p[i]; c3 <- pairs_df$cut_site_3p[i]
    cut_s <- min(c5, c3); cut_e <- max(c5, c3)
    del_size <- abs(c3 - c5)

    use_strat_b <- del_size > max_wt_amplicon

    if (!use_strat_b) {
      # ---- Strategy A: Flanking Pair ----
      strategies[i] <- "Flanking"

      tmpl_start <- max(1, cut_s - padding)
      tmpl_end   <- cut_e + padding

      seq_dna <- safe_get_seq(genome, names=chrom, start=tmpl_start, end=tmpl_end)

      if (nchar(seq_dna) > 0) {
        rel_start <- cut_s - tmpl_start
        rel_len   <- cut_e - cut_s

        req <- list(
          sequence_template = seq_dna,
          target_start = rel_start,
          target_len = rel_len,
          tm_opt = target_tm,
          prod_min = rel_len + 80,
          prod_max = rel_len + 1000
        )

        batch_requests[[length(batch_requests) + 1]] <- req
        pending_assignments[[length(pending_assignments) + 1]] <- list(
          row = i, type = "ext_stratA", batch_idx = length(batch_requests), del_size = del_size
        )
      }

    } else {
      # ---- Strategy B: Dual Pair ----
      strategies[i] <- "Dual_Pair"
      size_wt[i] <- "Too Large"

      # -- B1. External Pair (Mutant Detection) --
      up_seq <- safe_get_seq(genome, names=chrom, start=cut_s-400, end=cut_s)
      dn_seq <- safe_get_seq(genome, names=chrom, start=cut_e, end=cut_e+400)

      if (nchar(up_seq) > 0 && nchar(dn_seq) > 0) {
        mut_templ <- paste0(up_seq, dn_seq)

        req_ext <- list(
          sequence_template = mut_templ,
          target_start = 400,
          target_len = 0,
          tm_opt = target_tm,
          prod_min = 150,
          prod_max = 600
        )

        batch_requests[[length(batch_requests) + 1]] <- req_ext
        pending_assignments[[length(pending_assignments) + 1]] <- list(
          row = i, type = "ext_stratB", batch_idx = length(batch_requests)
        )
      }

      # -- B2. Internal Pair (WT Detection) --
      e5_rank <- pairs_df$exon_5p[i]; e3_rank <- pairs_df$exon_3p[i]
      cache_key <- paste(e5_rank, e3_rank, sep="_")

      if (!is.null(internal_cache[[cache_key]])) {
        # Hit!
        pending_assignments[[length(pending_assignments) + 1]] <- list(
          row = i, type = "int", batch_idx = internal_cache[[cache_key]]
        )
      } else {
        # Miss!
        del_exons <- exons_df[exons_df$rank > e5_rank & exons_df$rank < e3_rank, ]
        target_s <- NA; target_e <- NA

        # Priority 1: Largest deleted exon
        if (nrow(del_exons) > 0) {
          # FIX: Use 'end' and 'start' (GRanges defaults) instead of 'exon_chrom_end/start'
          # which are stripped when converting GRanges -> data.frame
          best_ex <- del_exons[which.max(del_exons$end - del_exons$start), ]
          w <- best_ex$end - best_ex$start

          # Ensure w is not empty/NA before comparison
          if (length(w) > 0 && !is.na(w) && w > 80) {
            target_s <- best_ex$start
            target_e <- best_ex$end
          }
        }

        # Priority 2: Midpoint of deletion
        if (is.na(target_s)) {
          mid <- floor((cut_s + cut_e) / 2)
          target_s <- mid - 200
          target_e <- mid + 200
        }

        int_seq <- safe_get_seq(genome, names=chrom, start=target_s, end=target_e)

        if (nchar(int_seq) > 100) {
          mid_idx <- floor(nchar(int_seq) / 2)

          req_int <- list(
            sequence_template = int_seq,
            target_start = mid_idx,
            target_len = 0,
            tm_opt = target_tm,
            prod_min = 80,
            prod_max = min(400, nchar(int_seq))
          )

          batch_requests[[length(batch_requests) + 1]] <- req_int
          new_idx <- length(batch_requests)
          internal_cache[[cache_key]] <- new_idx

          pending_assignments[[length(pending_assignments) + 1]] <- list(
            row = i, type = "int", batch_idx = new_idx
          )
        }
      }
    }
  }
  close(pb)

  # --- PHASE 2: Execute Batch in Python ---
  if (length(batch_requests) > 0) {
    message("Running Batch Primer3 (", length(batch_requests), " designs)...")
    batch_results <- tryCatch({
      run_primer3_batch(batch_requests)
    }, error = function(e) {
      warning("Python batch execution failed: ", e$message)
      return(NULL)
    })
  } else {
    batch_results <- NULL
  }

  # --- PHASE 3: Process Results ---
  if (!is.null(batch_results)) {
    message("Mapping results to dataframe...")

    for (assign in pending_assignments) {
      if (assign$batch_idx > length(batch_results)) next

      res <- batch_results[[assign$batch_idx]]
      r_idx <- assign$row

      if (!is.null(res)) {
        if (assign$type == "ext_stratA") {
          ext_fwd[r_idx] <- res$fwd_seq
          ext_rev[r_idx] <- res$rev_seq
          size_wt[r_idx] <- as.character(res$prod_size)
          size_mut[r_idx] <- res$prod_size - assign$del_size
        }
        else if (assign$type == "ext_stratB") {
          ext_fwd[r_idx] <- res$fwd_seq
          ext_rev[r_idx] <- res$rev_seq
          size_mut[r_idx] <- res$prod_size
        }
        else if (assign$type == "int") {
          int_fwd[r_idx] <- res$fwd_seq
          int_rev[r_idx] <- res$rev_seq
          size_int[r_idx] <- res$prod_size
        }
      }
    }
  }

  pairs_df$priming_strategy <- strategies
  pairs_df$primer_ext_fwd   <- ext_fwd
  pairs_df$primer_ext_rev   <- ext_rev
  pairs_df$primer_int_fwd   <- int_fwd
  pairs_df$primer_int_rev   <- int_rev
  pairs_df$exp_wt_size      <- size_wt
  pairs_df$exp_mut_size     <- size_mut
  pairs_df$exp_int_size     <- size_int

  n_succ <- sum(!is.na(pairs_df$primer_ext_fwd))
  message("Designed primers for ", n_succ, "/", n_total, " pairs.")

  return(pairs_df)
}
