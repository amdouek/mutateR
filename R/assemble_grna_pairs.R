#' Assemble and annotate gRNA pairs for exon‑flanking deletions
#'
#' Merges gRNAs for phase-compatible exon pairs.
#' Automatically detects scoring method to apply appropriate cutoffs:
#' - DeepCpf1 / DeepSpCas9: > 50 (linear regression on indel frequency)
#' - RuleSet1 / Azimuth / enPAM+GB: > 0.5 (probability-like scores)
#' - RuleSet3: > 0.1 (z-scores)
#'
#' @param grna_gr GRanges returned by \link{filter_valid_grnas}.
#' @param exon_gr GRanges from \link{get_exon_structures}(output="GRanges").
#' @param transcript_id Character. Ensembl transcript ID.
#' @param species Character. Ensembl species short code.
#' @param score_cutoff Numeric. Optional. If NULL, auto-selects based on `scoring_method` metadata.
#'
#' @return A \code{data.frame} of candidate exon‑flanking gRNA pairs.
#' @export
assemble_grna_pairs <- function(grna_gr,
                                exon_gr,
                                transcript_id,
                                species,
                                score_cutoff = NULL) {

  stopifnot(inherits(grna_gr, "GRanges"))
  message("Assembling gRNA pairs for exon‑flanking deletions...")

  # ---- 1. Determine score cutoff ------------------------------------
  if (is.null(score_cutoff)) {
    # Check metadata for scoring method
    method <- if ("scoring_method" %in% names(mcols(grna_gr))) unique(mcols(grna_gr)$scoring_method)[1] else "ruleset1"

    # Define model categories by output type
    # Regression models output raw indel percentages (0-100 scale)
    regression_models <- c("deepcpf1", "deepspcas9")

    # Z-score models output standardised activity scores (centered ~0)
    zscore_models <- c("ruleset3")

    # Probability-like models output normalised scores (~0-1 scale)
    # Note: enpamgb can slightly exceed 1.0 (or be negative) due to unconstrained regression
    probability_models <- c("ruleset1", "azimuth", "deephf", "enpamgb")

    if (!is.na(method) && tolower(method) %in% regression_models) {
      score_cutoff <- 50
      message("Detected linear regression scores (", method, "). Using default cutoff: ", score_cutoff)
    } else if (!is.na(method) && tolower(method) %in% zscore_models) {
      score_cutoff <- 0.1
      message("Detected z-scored activity model (", method, "). Using default cutoff: ", score_cutoff)
    } else {
      score_cutoff <- 0.5
      message("Detected probability-like scores (", method, "). Using default cutoff: ", score_cutoff)
    }
  }

  # ---- 2. Flatten on-target scores safely ---------------------------
  raw_scores <- mcols(grna_gr)$ontarget_score
  if (is.numeric(raw_scores)) {
    numeric_scores <- raw_scores
  } else {
    numeric_scores <- suppressWarnings(as.numeric(unlist(raw_scores)))
  }
  numeric_scores[is.nan(numeric_scores)] <- NA_real_

  # ---- 3. Intragenic mode (≤2 exons) --------------------------------
  n_exons <- length(exon_gr)
  if (n_exons <= 2) {
    message("Single-exon/two-exon gene detected: constructing intragenic deletion pairs.")
    site_df <- as.data.frame(grna_gr)
    if (!"cut_site" %in% names(site_df)) site_df$cut_site <- if("cut_site" %in% names(mcols(grna_gr))) mcols(grna_gr)$cut_site else site_df$start
    if (!"protospacer_sequence" %in% names(site_df)) site_df$protospacer_sequence <- substring(as.character(site_df$sequence_context), 5, nchar(as.character(site_df$sequence_context)) - 6)

    site_df$exon_rank <- if("exon_rank" %in% names(site_df)) site_df$exon_rank else 1L
    site_df$score <- numeric_scores

    if (nrow(site_df) < 2) return(NULL)

    comb_idx <- utils::combn(seq_len(nrow(site_df)), 2)
    intragenic <- data.frame(
      protospacer_sequence_5p = site_df$protospacer_sequence[comb_idx[1,]],
      protospacer_sequence_3p = site_df$protospacer_sequence[comb_idx[2,]],
      ontarget_score_5p       = site_df$score[comb_idx[1,]],
      ontarget_score_3p       = site_df$score[comb_idx[2,]],
      exon_5p                 = site_df$exon_rank[comb_idx[1,]],
      exon_3p                 = site_df$exon_rank[comb_idx[2,]],
      cut_site_5p             = site_df$cut_site[comb_idx[1,]],
      cut_site_3p             = site_df$cut_site[comb_idx[2,]],
      seqnames_5p             = site_df$seqnames[comb_idx[1,]],
      genomic_deletion_size   = abs(site_df$cut_site[comb_idx[2,]] - site_df$cut_site[comb_idx[1,]]),
      transcript_id           = transcript_id,
      stringsAsFactors = FALSE
    )

    intragenic$recommended <- with(intragenic,
                                   ontarget_score_5p >= score_cutoff &
                                     ontarget_score_3p >= score_cutoff)

    return(list(pairs = intragenic[intragenic$recommended == TRUE, ], intragenic_mode = TRUE))
  }

  # ---- 4. Multi-exon logic ------------------------------------------
  exon_meta <- as.data.frame(mcols(exon_gr))
  exon_meta$rank <- seq_len(nrow(exon_meta))
  comp_pairs <- check_exon_phase(exon_meta, include_contiguous = FALSE)
  comp_pairs <- subset(comp_pairs, compatible == TRUE)

  if (nrow(comp_pairs) == 0) return(NULL)

  # Check frame
  fs_list <- lapply(seq_len(nrow(comp_pairs)), function(i)
    with(comp_pairs[i, ], check_frameshift_ptc(exon_meta, exon_5p, exon_3p)))
  fs_df <- do.call(rbind, lapply(fs_list, as.data.frame))
  pair_info <- cbind(comp_pairs, fs_df)
  pair_info$target_5p <- pair_info$exon_5p + 1
  pair_info$target_3p <- pair_info$exon_3p - 1
  pair_info <- subset(pair_info, (compatible == TRUE | terminal_exon_case == TRUE) &
                        target_5p > 0 & target_3p <= nrow(exon_meta))

  if (nrow(pair_info) == 0) return(NULL)

  # Prepare dataframe
  grna_df <- as.data.frame(grna_gr)
  grna_df$ontarget_score <- numeric_scores
  if (!"cut_site" %in% names(grna_df)) grna_df$cut_site <- if("cut_site" %in% names(mcols(grna_gr))) mcols(grna_gr)$cut_site else grna_df$start
  if (!"exon_rank" %in% names(grna_df)) {
    hits <- GenomicRanges::findOverlaps(grna_gr, exon_gr)
    grna_df$exon_rank <- NA_integer_
    grna_df$exon_rank[queryHits(hits)] <- exon_meta$rank[subjectHits(hits)]
  }

  domain_df <- tryCatch(map_protein_domains(transcript_id, species), error = \(e) NULL)

  # Build pairs
  results <- list()
  for (i in seq_len(nrow(pair_info))) {
    e5 <- pair_info$target_5p[i]; e3 <- pair_info$target_3p[i]
    g5 <- subset(grna_df, exon_rank == e5)
    g3 <- subset(grna_df, exon_rank == e3)

    if (nrow(g5) > 0 && nrow(g3) > 0) {
      comb <- merge(g5, g3, by = NULL, suffixes = c("_5p", "_3p"))
      comb$upstream_pair <- pair_info$exon_5p[i]; comb$downstream_pair <- pair_info$exon_3p[i]
      comb$exon_5p <- e5; comb$exon_3p <- e3
      comb$compatible <- pair_info$compatible[i]
      comb$frameshift <- ifelse(pair_info$terminal_exon_case[i], TRUE, pair_info$frameshift[i])
      comb$ptc_flag   <- ifelse(pair_info$terminal_exon_case[i], FALSE, pair_info$ptc_flag[i])
      comb$terminal_exon_case <- pair_info$terminal_exon_case[i]
      comb$genomic_deletion_size <- abs(comb$cut_site_3p - comb$cut_site_5p)
      comb$transcript_deletion_size <- as.integer(pair_info$deleted_length[i])

      dom_overlap <- if (!is.null(domain_df)) subset(domain_df, exon_rank %in% seq(e5, e3)) else NULL
      comb$domains <- if (!is.null(dom_overlap) && nrow(dom_overlap) > 0) paste(unique(dom_overlap$domain_desc), collapse = "; ") else NA_character_

      results[[length(results) + 1]] <- comb
    }
  }

  if (!length(results)) return(NULL)
  out <- do.call(rbind, results)

  # Technical filter: Minimal deletion size & potential for steric hindrance
  # Deletions < 50 bp are hard to genotype by agarose gel
  # Cas effector footprint causing steric hindrance if two RNPs are too close to each other.
  out <- out[out$genomic_deletion_size > 50, ]

  if (nrow(out) > 0) {
    # Biological filter: Residual Exon Mass (REM)
    # For single-exon targets, ensure the exon is effectively destroyed/skipped.

    # Extract exon widths using the rank metadata
    ex_widths <- setNames(width(exon_gr), mcols(exon_gr)$rank)

    # Identify single-exon targets (Start Exon == End Exon)
    is_single_ex <- out$exon_5p == out$exon_3p

    if (any(is_single_ex)) {
      # Get width of the targeted exon
      target_w <- ex_widths[as.character(out$exon_5p)]
      del_s    <- out$genomic_deletion_size

      # Calculate biological metrics (residual exon length after deletion)
      residual_len <- target_w - del_s
      destruction_ratio <- del_s / target_w

      # Criteria for keeping single-exon targets:
      # A. Residual fragment is too small for spliceosome recognition (< 50 bp)
      # B. The deletion destroys the vast majority of the exon (> 70%)
      biologically_valid <- (residual_len < 50) | (destruction_ratio > 0.70)

      # Keep if: multi-exon deletion OR passes biological check
      to_keep <- (!is_single_ex) | biologically_valid

      out <- out[to_keep, ]
    }
  }

  if (nrow(out) == 0) return(NULL)

  # 'Recommended' logic
  out$recommended <- with(out,
                          !is.na(ontarget_score_5p) & ontarget_score_5p >= score_cutoff &
                            !is.na(ontarget_score_3p) & ontarget_score_3p >= score_cutoff)

  # Keep columns
  keep_cols <- c("upstream_pair","downstream_pair","exon_5p","exon_3p",
                 "compatible","frameshift","ptc_flag","terminal_exon_case",
                 "genomic_deletion_size", "transcript_deletion_size",
                 "protospacer_sequence_5p","pam_sequence_5p", "ontarget_score_5p",
                 "protospacer_sequence_3p","pam_sequence_3p", "ontarget_score_3p",
                 "domains","recommended",
                 "seqnames_5p", "start_5p", "end_5p", "cut_site_5p",
                 "seqnames_3p", "start_3p", "end_3p", "cut_site_3p")
  out <- out[, intersect(names(out), keep_cols), drop = FALSE]

  message("Generated ", nrow(out), " candidate exon‑flanking gRNA pairs.")
  return(out)
}
