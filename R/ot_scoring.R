#' @title Off-target scoring orchestrator and aggregation
#'
#' @description Backend-agnostic off-target scoring pipeline. The main
#' \code{score_offtargets()} function dispatches to Bowtie (fast, mismatch-only)
#' or CRISPRitz (slower, supports bulges) for genome-wide off-target search,
#' then scores and aggregates hits via the shared \code{score_and_aggregate()}.
#'
#' @keywords internal
#' @name ot_scoring
NULL


#' @title Off-target scoring for gRNAs
#'
#' @description Performs genome-wide off-target search using one of three backends,
#' scores each hit using the CFD method, and aggregates to a per-gRNA MIT-style
#' specificity score. All unique protospacers are searched once, then results are
#' mapped back to individual gRNAs with per-gRNA on-target exclusion.
#'
#' @param grna_gr GRanges. Returned by \code{\link{score_grnas}} (post on-target scoring).
#'        Must contain metadata columns: \code{protospacer_sequence}, \code{pam_sequence}, \code{cut_site}.
#' @param genome BSgenome object for the relevant species.
#' @param nuclease Character. One of "Cas9", "Cas12a", or "enCas12a".
#' @param exon_gr GRanges. Exon coordinates from \code{\link{get_exon_structures}}.
#'        Used for annotating off-target hits as exonic/non-exonic.
#' @param backend Character. Off-target search backend:
#'   \describe{
#'     \item{"bowtie"}{(Default) Fast mismatch-only search via crisprBowtie/Rbowtie.
#'           Typically completes in <1 min for ~400 guides at 3 mismatches.
#'           Requires a one-time Bowtie index build (~10–15 min for human genome).
#'           No bulge detection.}
#'     \item{"hybrid"}{Bowtie mismatch search on all guides (<1 min).
#'           Designed for use with \code{run_mutateR()}, which subsequently runs
#'           CRISPRitz bulge refinement on recommended-pair gRNAs only via
#'           \code{\link{refine_with_bulges}}. Falls back to Bowtie-only if
#'           CRISPRitz is unavailable.}
#'     \item{"crispritz"}{When called from \code{run_mutateR()}, runs a
#'           single-pass indexed search (bMax capped at 1) on all guides,
#'           capturing both mismatch and bulge off-targets. No Bowtie
#'           dependency. When called directly, searches with the exact
#'           parameters provided. Requires CRISPRitz in the mutateR conda
#'           environment (or WSL on Windows).}
#'   }
#' @param max_mismatches Integer. Maximum allowed mismatches (default 3).
#'        Clamped to 3 for Bowtie backend (hard limit). CRISPRitz supports up to 6.
#' @param max_dna_bulges Integer. Maximum DNA bulge size (default 2).
#'        Only used by CRISPRitz/hybrid backends.
#' @param max_rna_bulges Integer. Maximum RNA bulge size (default 2).
#'        Only used by CRISPRitz/hybrid backends.
#' @param bulge_penalty Numeric. Multiplicative penalty per bulge nucleotide applied
#'        to CFD scores (default 0.2). Only affects CRISPRitz results containing bulges.
#' @param canonical_only Logical (default TRUE). For Bowtie backend, if TRUE return
#'        only hits adjacent to canonical PAM sequences. The MIT specificity score was
#'        calibrated against canonical PAM off-targets (Cas9 NGG). Set to FALSE for
#'        comprehensive scanning of non-canonical PAMs; note that this drastically
#'        increases hit volume and reduces specificity scores. Ignored by CRISPRitz
#'        backend (PAM is specified in the PAM file).
#' @param detail_level Character. Controls how much off-target hit detail is
#'   retained in the output (default \code{"compact"}). Per-gRNA specificity
#'   scores are always computed from the \strong{complete} hit set regardless
#'   of this setting; only the attached detail table is affected.
#'   \describe{
#'     \item{"compact"}{(Default) Retains hits with CFD > 0.01, capped at
#'           200 per gRNA. Suitable for inspection and reporting.}
#'     \item{"full"}{Retains all hits (may exceed 1 GB for large genomes).
#'           Use for custom re-analysis of the full off-target landscape.}
#'     \item{"file"}{Writes full details to a compressed RDS file in
#'           \code{cache_dir}; attaches compact version in memory. File
#'           path stored as \code{attr(result, "offtarget_details_path")}.}
#'     \item{"none"}{Discards all hit details. Smallest memory footprint;
#'           per-gRNA scores are still computed and attached.}
#'   }
#' @param threads Integer. Number of threads for CRISPRitz (default 4). Ignored by Bowtie.
#' @param timeout Integer. Maximum seconds for CRISPRitz execution (default 600).
#'        Ignored by Bowtie.
#' @param cache_dir Character or NULL. Directory for cached genome data (FASTA + indices).
#'        Defaults to \code{tools::R_user_dir("mutateR", "cache")}.
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return The input GRanges with additional metadata columns:
#'   \describe{
#'     \item{specificity_score}{Numeric (0-100). MIT-style aggregate specificity.
#'           100 = no off-targets found; lower = more off-target risk.}
#'     \item{n_offtargets}{Integer. Number of off-target hits with CFD > 0.01.}
#'     \item{has_exonic_offtarget}{Logical. Whether any exonic off-target has CFD > 0.1.}
#'   }
#'   The off-target hit table (filtered per \code{detail_level}) is attached as
#'   \code{attr(result, "offtarget_details")}. When \code{detail_level = "file"},
#'   the full table filepath is also stored as
#'   \code{attr(result, "offtarget_details_path")}.
#'   The backend used is recorded as \code{attr(result, "ot_backend")}.
#'
#' @seealso \code{\link{score_grnas}} for on-target scoring,
#'          \code{\link{ensure_bowtie_index}} for Bowtie index building,
#'          \code{\link{install_mutater_env}} for CRISPRitz installation.
#'
#' @export
score_offtargets <- function(grna_gr,
                             genome,
                             nuclease = c("Cas9", "Cas12a", "enCas12a"),
                             exon_gr,
                             backend = c("bowtie", "hybrid", "crispritz"),
                             max_mismatches = 3L,
                             max_dna_bulges = 2L,
                             max_rna_bulges = 2L,
                             bulge_penalty = 0.2,
                             canonical_only = TRUE,
                             detail_level = c("compact", "full", "file", "none"),
                             threads = 4L,
                             timeout = 600L,
                             cache_dir = NULL,
                             quiet = FALSE) {

  nuclease <- match.arg(nuclease)
  backend  <- match.arg(backend)
  detail_level <- match.arg(detail_level)

  # ---- 1. Input validation (shared) ----
  if (!inherits(grna_gr, "GRanges")) {
    stop("'grna_gr' must be a GRanges object.")
  }

  required_mcols <- c("protospacer_sequence", "pam_sequence", "cut_site")
  missing_mcols <- setdiff(required_mcols, names(mcols(grna_gr)))
  if (length(missing_mcols) > 0) {
    stop("GRanges is missing required metadata columns: ",
         paste(missing_mcols, collapse = ", "),
         "\nEnsure grna_gr was produced by find_cas9_sites() or find_cas12a_sites().")
  }

  if (!inherits(genome, "BSgenome")) {
    stop("'genome' must be a BSgenome object.")
  }

  if (!inherits(exon_gr, "GRanges")) {
    stop("'exon_gr' must be a GRanges object.")
  }

  if (!is.numeric(max_mismatches) || max_mismatches < 0 || max_mismatches > 6) {
    stop("'max_mismatches' must be an integer between 0 and 6.")
  }
  if (!is.numeric(bulge_penalty) || bulge_penalty <= 0 || bulge_penalty >= 1) {
    stop("'bulge_penalty' must be a numeric value in (0, 1).")
  }

  # Backend-specific validation
  if (backend == "crispritz") {
    if (!is.numeric(max_dna_bulges) || max_dna_bulges < 0 || max_dna_bulges > 3) {
      stop("'max_dna_bulges' must be an integer between 0 and 3.")
    }
    if (!is.numeric(max_rna_bulges) || max_rna_bulges < 0 || max_rna_bulges > 3) {
      stop("'max_rna_bulges' must be an integer between 0 and 3.")
    }
    if (!check_mutater_env()) {
      stop("mutateR Python environment is not available.\n",
           "Please run install_mutater_env() and restart R before using the CRISPRitz backend.\n",
           "Alternatively, use backend = 'bowtie' for mismatch-only off-target scoring.")
    }
  }

  n_grnas <- length(grna_gr)
  if (!quiet) message("Off-target analysis for ", n_grnas, " gRNAs (backend: ", backend, ")...")

  # ---- 2. Setup (shared) ----
  if (is.null(cache_dir)) {
    cache_dir <- tools::R_user_dir("mutateR", "cache")
  }
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- 3. Deduplicate protospacers (shared) ----
  dedup <- deduplicate_protospacers(grna_gr)
  guide_map      <- dedup$guide_map
  unique_spacers <- dedup$unique_spacers

  if (!quiet) {
    message("Unique protospacers: ", length(unique_spacers),
            " (from ", n_grnas, " total gRNAs; ",
            n_grnas - length(unique_spacers), " duplicates removed for search)")
  }

  # ---- 4. Backend dispatch ----
  t_start <- Sys.time()

  if (backend %in% c("bowtie", "hybrid")) {

    hit_df <- run_bowtie_backend(
      unique_spacers = unique_spacers,
      guide_map      = guide_map,
      genome         = genome,
      nuclease       = nuclease,
      max_mismatches = max_mismatches,
      canonical_only = canonical_only,
      cache_dir      = cache_dir,
      quiet          = quiet
    )

  } else {
    # backend == "crispritz"
    hit_df <- run_crispritz_backend(
      unique_spacers = unique_spacers,
      guide_map      = guide_map,
      genome         = genome,
      nuclease       = nuclease,
      max_mismatches = max_mismatches,
      max_dna_bulges = max_dna_bulges,
      max_rna_bulges = max_rna_bulges,
      threads        = threads,
      timeout        = timeout,
      cache_dir      = cache_dir,
      quiet          = quiet
    )
  }

  t_elapsed <- difftime(Sys.time(), t_start, units = "secs")
  if (!quiet) message("Off-target search completed in ",
                      round(as.numeric(t_elapsed), 1), " seconds.")

  if (!quiet) message("Found ", nrow(hit_df), " total hits ",
                      "(including on-target loci; these will be excluded during aggregation).")

  # ---- 5a. Score and aggregate (shared) ----
  if (!quiet) message("Scoring off-target hits and computing specificity...")
  agg <- score_and_aggregate(
    hit_df        = hit_df,
    grna_gr       = grna_gr,
    exon_gr       = exon_gr,
    nuclease      = nuclease,
    bulge_penalty = bulge_penalty
  )

  # ---- 5a-i. Mark unsearched gRNAs (partial batch failure) ----
  # When the batched CRISPRitz backend suffers partial failures (OOM, WSL
  # degradation), some spacers are never searched. score_and_aggregate()
  # assigns them specificity=100 (no hits = no off-target burden), which is
  # an artifact — absence of evidence ≠ evidence of absence.
  # Replace with NA so downstream consumers can distinguish "genuinely
  # unique" from "unsearched".
  all_spacers <- toupper(as.character(mcols(grna_gr)$protospacer_sequence))

  searched_spacers <- attr(hit_df, "searched_spacers")
  if (!is.null(searched_spacers)) {
    unsearched   <- !all_spacers %in% toupper(searched_spacers)
    n_unsearched <- sum(unsearched)

    if (n_unsearched > 0) {
      agg$per_grna$specificity_score[unsearched]    <- NA_real_
      agg$per_grna$n_offtargets[unsearched]         <- NA_integer_
      agg$per_grna$has_exonic_offtarget[unsearched]  <- NA

      if (!quiet) {
        n_searched <- sum(!unsearched)
        message("Note: ", n_unsearched, " of ", n_grnas, " gRNA(s) could not be ",
                "searched (batch failures). Their specificity is set to NA.")
        message("  Successfully searched: ", n_searched, " gRNAs (",
                length(unique(searched_spacers)), " unique spacers)")
      }
    }
  }

  # ---- 5a-ii. Mark quarantined gRNAs (repeat-rich, known high off-target risk) ----
  quarantined_spacers <- attr(hit_df, "quarantined_spacers")
  if (!is.null(quarantined_spacers) && length(quarantined_spacers) > 0) {
    is_quarantined <- all_spacers %in% toupper(quarantined_spacers)
    n_quarantined  <- sum(is_quarantined)

    if (n_quarantined > 0) {
      # Specificity = 0: known high off-target risk (not unknown like NA)
      agg$per_grna$specificity_score[is_quarantined]    <- 0
      # n_offtargets and exonic status are genuinely unknown
      agg$per_grna$n_offtargets[is_quarantined]         <- NA_integer_
      agg$per_grna$has_exonic_offtarget[is_quarantined] <- NA

      if (!quiet) {
        message("Note: ", n_quarantined, " gRNA(s) quarantined as repeat-rich ",
                "(CRISPRitz OOM). Specificity set to 0.")
      }
    }
  }

  # ---- 5b. Apply detail level compaction ----
  ot_details      <- NULL
  ot_details_path <- NULL

  if (detail_level == "full") {
    ot_details <- agg$hit_details

  } else if (detail_level == "compact") {
    ot_details <- compact_hit_details(agg$hit_details)
    if (!quiet) {
      n_full    <- nrow(agg$hit_details)
      n_compact <- if (is.null(ot_details)) 0L else nrow(ot_details)
      message("Off-target details compacted: ", n_full, " -> ", n_compact, " hits retained")
    }

  } else if (detail_level == "file") {
    # Write full details to compressed RDS
    ot_details_path <- file.path(
      cache_dir,
      paste0("offtarget_details_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    )
    saveRDS(agg$hit_details, ot_details_path, compress = "gzip")
    if (!quiet) message("Full off-target details written to: ", ot_details_path,
                        " (", round(file.size(ot_details_path) / 1e6, 1), " MB)")
    ot_details <- compact_hit_details(agg$hit_details)

  }
  # detail_level == "none" -- ot_details stays NULL

  # ---- 6. Attach results to GRanges (shared) ----
  mcols(grna_gr)$specificity_score    <- agg$per_grna$specificity_score
  mcols(grna_gr)$n_offtargets         <- agg$per_grna$n_offtargets
  mcols(grna_gr)$has_exonic_offtarget <- agg$per_grna$has_exonic_offtarget

  attr(grna_gr, "offtarget_details")      <- ot_details
  attr(grna_gr, "offtarget_details_path") <- ot_details_path
  attr(grna_gr, "ot_backend")             <- backend
  attr(grna_gr, "ot_detail_level")        <- detail_level

  # ---- 7. Summary reporting (shared) ----
  if (!quiet) {
    median_spec <- median(agg$per_grna$specificity_score, na.rm = TRUE)
    n_high      <- sum(agg$per_grna$specificity_score >= 50, na.rm = TRUE)
    n_low       <- sum(agg$per_grna$specificity_score < 50, na.rm = TRUE)
    n_exonic    <- sum(agg$per_grna$has_exonic_offtarget, na.rm = TRUE)

    message("\n--- Off-target scoring summary ---")
    message("  Backend:               ", backend,
            if (backend == "hybrid") " (mismatch pass; bulge refinement pending)" else "")
    message("  Total gRNAs scored:    ", n_grnas)
    message("  Median specificity:    ", round(median_spec, 1))
    message("  High specificity (>=50): ", n_high, " (", round(100 * n_high / n_grnas, 1), "%)")
    message("  Low specificity (<50):   ", n_low, " (", round(100 * n_low / n_grnas, 1), "%)")
    if (n_exonic > 0) {
      message("  WARNING: ", n_exonic, " gRNA(s) have exonic off-target hits (CFD > 0.1)")
    }
    message("---------------------------------\n")
  }

  return(grna_gr)
}


# ----- Backend dispatch helpers (internal) -----


#' Run the Bowtie off-target search pipeline
#'
#' Orchestrates index caching, search, and normalisation for the Bowtie backend.
#' Called by \code{score_offtargets()} when \code{backend} is \code{"bowtie"} or
#' \code{"hybrid"}.
#'
#' @param unique_spacers Character vector. Deduplicated protospacer sequences.
#' @param guide_map Named list. Protospacer → GRanges index mapping.
#' @param genome BSgenome object.
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a".
#' @param max_mismatches Integer.
#' @param canonical_only Logical.
#' @param cache_dir Character. Cache directory.
#' @param quiet Logical.
#'
#' @return Standardised hit_df data.frame (same schema as \code{parse_crispritz_output}).
#' @keywords internal
run_bowtie_backend <- function(unique_spacers,
                               guide_map,
                               genome,
                               nuclease,
                               max_mismatches,
                               canonical_only = TRUE,
                               cache_dir,
                               quiet) {

  # ---- 1. Ensure Bowtie index (cascades to FASTA cache) ----
  if (!quiet) message("Checking Bowtie index cache...")
  bowtie_index <- ensure_bowtie_index(genome, cache_dir, quiet = quiet)

  # ---- 2. Run search ----
  if (!quiet) message("Launching Bowtie genome-wide search...")

  bowtie_raw <- run_bowtie_search(
    unique_spacers = unique_spacers,
    genome         = genome,
    nuclease       = nuclease,
    bowtie_index   = bowtie_index,
    max_mismatches = max_mismatches,
    canonical_only = canonical_only,
    quiet          = quiet
  )

  # ---- 3. Normalise to standard hit_df schema ----
  if (!quiet) message("Normalising Bowtie hits...")
  hit_df <- normalise_bowtie_hits(bowtie_raw, guide_map, nuclease)

  return(hit_df)
}


#' Run the CRISPRitz off-target search pipeline
#'
#' Orchestrates FASTA caching, input preparation, CLI execution, and output
#' parsing for the CRISPRitz backend. Called by \code{score_offtargets()}
#' when \code{backend} is \code{"crispritz"}.
#'
#' @param unique_spacers Character vector. Deduplicated protospacer sequences.
#' @param guide_map Named list. Protospacer → GRanges index mapping.
#' @param genome BSgenome object.
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a".
#' @param max_mismatches Integer.
#' @param max_dna_bulges Integer.
#' @param max_rna_bulges Integer.
#' @param threads Integer.
#' @param timeout Integer.
#' @param cache_dir Character. Cache directory.
#' @param quiet Logical.
#'
#' @return Standardised hit_df data.frame.
#' @keywords internal
run_crispritz_backend <- function(unique_spacers,
                                  guide_map,
                                  genome,
                                  nuclease,
                                  max_mismatches,
                                  max_dna_bulges,
                                  max_rna_bulges,
                                  threads,
                                  timeout,
                                  cache_dir,
                                  quiet) {

  # ---- 1a. FASTA cache ----
  if (!quiet) message("Checking genome FASTA cache...")
  genome_dir <- ensure_genome_cache(genome, cache_dir, quiet = quiet)

  # ---- 1b. Ensure CRISPRitz genome index (bulge mode only) ----
  # Note: ensure_crispritz_index() uses its own default timeout (3600s),
  # independent of the search timeout. Genome indexing is a one-time
  # operation that typically takes 20-60 min for hg38; the search timeout
  # (default 600s) is far too short and would cause spurious failures.
  max_bulge  <- max(as.integer(max_dna_bulges), as.integer(max_rna_bulges))
  cz_index   <- ensure_crispritz_index(
    genome    = genome,
    nuclease  = nuclease,
    max_bulge = max_bulge,
    cache_dir = cache_dir,
    threads   = threads,
    quiet     = quiet
  )

  # ---- 1c. Cap threads for indexed bulge search ----
  # Empirical finding: CRISPRitz indexed search with >4 threads shows
  # severe performance regression (3x slower at 8 threads) and hit count
  # inflation, likely due to thread-unsafe TST traversal. Cap at 4.
  effective_threads <- threads
  if (cz_index$indexed && max_bulge > 0 && threads > 4L) {
    effective_threads <- 4L
    if (!quiet) message("Note: Capping threads to 4 for indexed bulge search ",
                        "(higher thread counts cause performance regression).")
  }

  # ---- 2. Route: batched indexed search vs single-call ----
  # CRISPRitz indexed bulge search pre-allocates per-guide data structures
  # and accumulates hits in memory. Large guide sets (>50) cause OOM-kill
  # in memory-constrained environments (e.g. WSL 14.7 GB). Split into
  # batches to keep per-invocation memory bounded.
  # Non-indexed (brute-force FASTA) mode does not have this problem —
  # it streams through chromosomes sequentially.
  use_batching <- cz_index$indexed && max_bulge > 0 &&
    length(unique_spacers) > 50L

  if (use_batching) {
    if (!quiet) message("Large guide set (", length(unique_spacers),
                        ") with indexed bulge search: using batched mode.")

    hit_df <- run_crispritz_batched(
      unique_spacers = unique_spacers,
      guide_map      = guide_map,
      genome_dir     = genome_dir,
      index_dir      = cz_index$index_dir,
      nuclease       = nuclease,
      max_mismatches = as.integer(max_mismatches),
      max_dna_bulges = as.integer(max_dna_bulges),
      max_rna_bulges = as.integer(max_rna_bulges),
      threads        = as.integer(effective_threads),
      timeout        = as.integer(timeout),
      batch_size     = 20L,
      quiet          = quiet
    )

    return(hit_df)
  }

  # ---- 3. Single-call path (small guide sets or non-indexed mode) ----
  tmpdir <- file.path(tempdir(), paste0("mutateR_ot_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

  # Preserve tmpdir for debugging; revert to on.exit(unlink(...)) when stable
  if (!quiet) {
    on.exit(message("CRISPRitz working directory: ", tmpdir), add = TRUE)
  }

  # ---- 4. Write inputs ----
  if (!quiet) message("Preparing CRISPRitz input files...")
  input_files <- write_crispritz_inputs(unique_spacers, nuclease, tmpdir)

  # Reconcile guide_map if input sanitisation removed sequences
  if (length(input_files$clean_seqs) < length(unique_spacers)) {
    removed_seqs <- setdiff(unique_spacers, input_files$clean_seqs)
    if (!quiet && length(removed_seqs) > 0) {
      message("Note: ", length(removed_seqs),
              " protospacer(s) removed during input sanitisation.")
    }
    guide_map <- guide_map[input_files$clean_seqs]
  }

  # ---- 5. Run search ----
  if (!quiet) {
    use_bulges  <- max_dna_bulges > 0 || max_rna_bulges > 0
    est_factor  <- if (use_bulges && cz_index$indexed) 0.15 else 0.5
    est_seconds <- length(unique_spacers) * max_mismatches * est_factor
    est_minutes <- round(est_seconds / 60, 1)
    message("Estimated CRISPRitz search time: ~", est_minutes, " minutes ",
            "(", length(input_files$clean_seqs), " guides, ",
            max_mismatches, "mm",
            if (use_bulges) paste0(", bDNA=", max_dna_bulges,
                                   ", bRNA=", max_rna_bulges) else " mismatch-only",
            ", ", effective_threads, " threads)")
    message("Launching CRISPRitz genome-wide search...")
  }

  result_file <- run_crispritz_search(
    genome_dir     = genome_dir,
    input_files    = input_files,
    index_dir      = if (cz_index$indexed) cz_index$index_dir else NULL,
    max_mismatches = as.integer(max_mismatches),
    max_dna_bulges = as.integer(max_dna_bulges),
    max_rna_bulges = as.integer(max_rna_bulges),
    threads        = as.integer(effective_threads),
    timeout        = as.integer(timeout),
    quiet          = quiet
  )

  # ---- 6. Parse output ----
  if (!quiet) message("Parsing CRISPRitz off-target hits...")
  hit_df <- parse_crispritz_output(
    result_file,
    guide_map,
    pam_length = input_files$pam_length,
    pam_side   = input_files$pam_side
  )

  return(hit_df)
}

#' Execute a single CRISPRitz write → search → parse cycle
#'
#' Encapsulates the core CRISPRitz workflow for a batch of guides:
#' write input files, run CLI search, parse output. Used by both the
#' single-call path in \code{\link{run_crispritz_backend}} and the
#' batching layer in \code{\link{run_crispritz_batched}}.
#'
#' @param spacers Character vector. Protospacer sequences for this chunk.
#' @param guide_map Named list. Subset of full guide_map for these spacers.
#' @param genome_dir Character. Path to cached FASTA directory.
#' @param index_dir Character or NULL. Path to CRISPRitz TST index directory.
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a".
#' @param max_mismatches Integer.
#' @param max_dna_bulges Integer.
#' @param max_rna_bulges Integer.
#' @param threads Integer.
#' @param timeout Integer. Seconds.
#' @param crispritz_info Named list or NULL. Pre-resolved CRISPRitz executable
#'   info, passed through to \code{\link{run_crispritz_search}}.
#' @param quiet Logical.
#'
#' @return data.frame with the same schema as \code{parse_crispritz_output()}.
#'   Returns a zero-row data.frame on empty output. Throws on search failure
#'   (caller should wrap in \code{tryCatch}).
#' @keywords internal
run_crispritz_chunk <- function(spacers, guide_map, genome_dir, index_dir,
                                nuclease, max_mismatches, max_dna_bulges,
                                max_rna_bulges, threads, timeout,
                                crispritz_info = NULL, quiet) {

  # Unique tmpdir per chunk to avoid file collisions in batched mode
  tmpdir <- file.path(tempdir(), paste0(
    "mutateR_ot_", format(Sys.time(), "%Y%m%d%H%M%S"),
    "_", sprintf("%04d", sample.int(9999, 1))
  ))
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  # ---- Write inputs ----
  input_files <- write_crispritz_inputs(spacers, nuclease, tmpdir)

  # Reconcile guide_map if sanitisation removed sequences
  if (length(input_files$clean_seqs) < length(spacers)) {
    guide_map <- guide_map[input_files$clean_seqs]
  }

  # ---- Search ----
  result_file <- run_crispritz_search(
    genome_dir     = genome_dir,
    input_files    = input_files,
    index_dir      = index_dir,
    max_mismatches = as.integer(max_mismatches),
    max_dna_bulges = as.integer(max_dna_bulges),
    max_rna_bulges = as.integer(max_rna_bulges),
    threads        = as.integer(threads),
    crispritz_info = crispritz_info,
    timeout        = as.integer(timeout),
    quiet          = quiet
  )

  # ---- Parse ----
  hit_df <- parse_crispritz_output(
    result_file,
    guide_map,
    pam_length = input_files$pam_length,
    pam_side   = input_files$pam_side
  )

  return(hit_df)
}

#' Probe individual guides from a failed CRISPRitz batch
#'
#' Replaces a failed batch with individual guide-by-guide probing to
#' identify and quarantine guides that crash CRISPRitz (repeat-rich OOM)
#' or produce excessive output (grey-zone guides matching millions of
#' genomic loci). Non-problematic guides are searched successfully and
#' their results saved.
#'
#' Each guide is run as a separate CRISPRitz invocation with a tight
#' per-guide timeout. Failure modes (OOM crash, output size exceeding
#' \code{mutateR.max_crispritz_output_bytes}, timeout, any other error)
#' all result in immediate quarantine with \code{specificity_score = 0}.
#'
#' More efficient than recursive bisection for this problem because:
#' (1) normal guides complete in ~49s each (predictable), (2) toxic
#' guides fail instantly (OOM or size guard), and (3) no cascading
#' WSL crashes from sub-batches containing hidden toxic guides.
#'
#' @param spacers Character vector. Protospacer sequences from the failed batch.
#' @param guide_map Named list. Subset of full guide_map for these spacers.
#' @param genome_dir Character. Path to cached FASTA directory.
#' @param index_dir Character. Path to CRISPRitz TST index directory.
#' @param nuclease Character. One of \code{"Cas9"}, \code{"Cas12a"}, \code{"enCas12a"}.
#' @param max_mismatches Integer. Maximum allowed mismatches.
#' @param max_dna_bulges Integer. Maximum DNA bulge size.
#' @param max_rna_bulges Integer. Maximum RNA bulge size.
#' @param threads Integer. Thread count for CRISPRitz.
#' @param timeout Integer. Caller's batch timeout (used for context only).
#'   Per-guide probe timeout is derived from the empirical scaling model:
#'   \code{(34 + 15) * 2.5 = 122s}, rounded to 120s.
#' @param crispritz_info Named list. Pre-resolved CRISPRitz executable info from
#'   \code{\link{locate_crispritz_bin}}, passed through to avoid re-probing
#'   degraded WSL.
#' @param batch_tmpdir Character. Directory for RDS result files.
#' @param batch_id Integer. Parent batch number, used in RDS file naming
#'   (pattern: \code{batch_{batch_id}_g{guide_index}.rds}).
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return A named list:
#'   \describe{
#'     \item{batch_files}{Character vector. Paths to RDS files containing
#'           successful single-guide results.}
#'     \item{searched_spacers}{Character vector. Spacers successfully searched.}
#'     \item{quarantined_spacers}{Character vector. Spacers that failed —
#'           confirmed problematic. Downstream consumers assign
#'           \code{specificity_score = 0}.}
#'     \item{total_hits}{Integer. Cumulative hit count across successful probes.}
#'   }
#'
#' @seealso \code{\link{run_crispritz_batched}} (sole caller),
#'   \code{\link{run_crispritz_chunk}} (atomic search unit),
#'   \code{\link{screen_repeat_rich}} for proactive pre-screening.
#'
#' @keywords internal
probe_failed_batch <- function(spacers, guide_map, genome_dir, index_dir,
                               nuclease, max_mismatches, max_dna_bulges,
                               max_rna_bulges, threads, timeout,
                               crispritz_info, batch_tmpdir, batch_id,
                               quiet = FALSE) {

  batch_files         <- character(0)
  searched_spacers    <- character(0)
  quarantined_spacers <- character(0)
  total_hits          <- 0L

  if (length(spacers) == 0) {
    return(list(batch_files = batch_files, searched_spacers = searched_spacers,
                quarantined_spacers = quarantined_spacers, total_hits = total_hits))
  }

  # Per-guide timeout: empirical model (34s fixed + 15s/guide) × 2.5 safety
  # Deliberately tight to catch grey-zone guides that don't OOM but produce
  # massive output. Normal guides finish in ~49s.
  probe_timeout <- 120L

  # Consecutive failure tracking: 5 in a row suggests system degradation
  # (not guide-specific), so quarantine the rest to avoid wasting time.
  max_consec_fail <- 5L
  consec_fail     <- 0L

  if (!quiet) message("    Probing ", length(spacers), " guides individually ",
                      "(timeout: ", probe_timeout, "s per guide)...")

  for (j in seq_along(spacers)) {
    spacer_j <- spacers[j]

    # ---- WSL health check with restart ----
    if (.Platform$OS.type == "windows") {
      wsl_ok <- tryCatch({
        res <- system2("wsl", args = c("bash", "-lc", shQuote("echo ok")),
                       stdout = TRUE, stderr = TRUE)
        exit <- attr(res, "status")
        if (is.null(exit)) exit <- 0L
        exit == 0L && any(grepl("ok", res))
      }, error = function(e) FALSE)

      if (!wsl_ok) {
        if (!quiet) message("      WSL unresponsive before guide ", j,
                            ". Restarting...")
        tryCatch(
          system2("wsl", args = "--shutdown", stdout = TRUE, stderr = TRUE),
          error = function(e) NULL
        )
        Sys.sleep(15)

        wsl_ok <- tryCatch({
          res <- system2("wsl", args = c("bash", "-lc", shQuote("echo ok")),
                         stdout = TRUE, stderr = TRUE)
          exit <- attr(res, "status")
          if (is.null(exit)) exit <- 0L
          exit == 0L && any(grepl("ok", res))
        }, error = function(e) FALSE)

        if (!wsl_ok) {
          n_remaining <- length(spacers) - j + 1L
          if (!quiet) message("      WSL failed to restart. Quarantining remaining ",
                              n_remaining, " guides.")
          quarantined_spacers <- c(quarantined_spacers, spacers[j:length(spacers)])
          break
        }
      }
    }

    # ---- Probe single guide ----
    probe_result <- tryCatch(
      run_crispritz_chunk(
        spacers        = spacer_j,
        guide_map      = guide_map[spacer_j],
        genome_dir     = genome_dir,
        index_dir      = index_dir,
        nuclease       = nuclease,
        max_mismatches = max_mismatches,
        max_dna_bulges = max_dna_bulges,
        max_rna_bulges = max_rna_bulges,
        threads        = threads,
        timeout        = probe_timeout,
        crispritz_info = crispritz_info,
        quiet          = TRUE
      ),
      error = function(e) {
        if (!quiet) message("      Guide ", j, " (", substr(spacer_j, 1, 12),
                            "...): FAILED — ", sub("\n.*", "", e$message))
        NULL
      }
    )

    # ---- Handle result ----
    if (!is.null(probe_result)) {
      n_hits <- nrow(probe_result)

      if (n_hits > 0) {
        rds_path <- file.path(batch_tmpdir,
                              paste0("batch_", batch_id, "_g", j, ".rds"))
        saveRDS(probe_result, rds_path)
        batch_files <- c(batch_files, rds_path)
      }

      searched_spacers <- c(searched_spacers, spacer_j)
      total_hits       <- total_hits + n_hits
      consec_fail      <- 0L

      if (!quiet) message("      Guide ", j, " (", substr(spacer_j, 1, 12),
                          "...): OK (", format(n_hits, big.mark = ","), " hits)")
    } else {
      quarantined_spacers <- c(quarantined_spacers, spacer_j)
      consec_fail         <- consec_fail + 1L

      if (consec_fail >= max_consec_fail) {
        n_remaining <- length(spacers) - j
        if (!quiet) message("      ", max_consec_fail, " consecutive failures. ",
                            "Quarantining remaining ", n_remaining, " guides.")
        if (j < length(spacers)) {
          quarantined_spacers <- c(quarantined_spacers, spacers[(j + 1L):length(spacers)])
        }
        break
      }
    }

    # ---- Inter-guide cleanup ----
    rm(probe_result)
    gc(full = TRUE, verbose = FALSE)
    if (j < length(spacers)) {
      Sys.sleep(1)
      reclaim_wsl_memory(quiet = TRUE)
    }
  }

  if (!quiet) {
    message("    Probe complete: ", length(searched_spacers), " searched, ",
            length(quarantined_spacers), " quarantined, ",
            format(total_hits, big.mark = ","), " total hits")
  }

  list(
    batch_files         = batch_files,
    searched_spacers    = searched_spacers,
    quarantined_spacers = quarantined_spacers,
    total_hits          = total_hits
  )
}

#' Batched CRISPRitz indexed search for large guide sets
#'
#' Splits guides into manageable batches to avoid OOM during indexed
#' bulge search. CRISPRitz pre-allocates per-guide data structures and
#' accumulates hits in memory before flushing; with many guides this
#' exceeds available RAM (empirically: 399 guides at bMax=1 triggers
#' OOM-kill in WSL at 14.7 GB).
#'
#' Each batch runs as a separate CRISPRitz invocation. Failed batches
#' are probed guide-by-guide via \code{\link{probe_failed_batch}} to
#' isolate and quarantine repeat-rich guides. Results are concatenated
#' across batches — no cross-batch deduplication is needed since each
#' batch searches disjoint guide sets.
#'
#' @param unique_spacers Character vector. Full set of unique protospacers.
#' @param guide_map Named list. Protospacer → GRanges index mapping.
#' @param genome_dir Character. Path to cached FASTA directory.
#' @param index_dir Character. Path to CRISPRitz TST index directory.
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a".
#' @param max_mismatches Integer.
#' @param max_dna_bulges Integer.
#' @param max_rna_bulges Integer.
#' @param threads Integer.
#' @param timeout Integer. Passed from caller; used only as context — actual
#'   per-batch timeouts are computed from the empirical scaling model.
#' @param batch_size Integer. Maximum guides per batch (default 20).
#'   Empirically tuned for CRISPRitz bMax=1 with 3 threads on hg38 in
#'   memory-constrained environments. Smaller batches keep per-invocation
#'   CRISPRitz memory bounded AND reduce R-side peak allocation
#'   during parsing (even with disk dedup, each batch produces a
#'   ~50-100 MB data.frame). With disk dedup active,
#'   per-batch overhead is ~1 min, so 20 batches vs 14 adds ~6 min.
#' @param quiet Logical.
#'
#' @return Concatenated hit_df from all batches (same schema as
#'   \code{parse_crispritz_output()}). Returns zero-row data.frame if
#'   all batches fail.
#' @keywords internal
run_crispritz_batched <- function(unique_spacers, guide_map, genome_dir,
                                  index_dir, nuclease, max_mismatches,
                                  max_dna_bulges, max_rna_bulges,
                                  threads, timeout,
                                  batch_size = 20L,
                                  quiet = FALSE) {

  n_spacers <- length(unique_spacers)
  n_batches <- ceiling(n_spacers / batch_size)

  # Per-batch timeout from empirical scaling model (bMax=1, 4 threads, hg38):
  #   time = 34s fixed + 15s/guide; ×2.5 safety margin.
  # Floor at 300s to accommodate system variability on small batches.
  batch_timeout <- as.integer(max(
    300L,
    ceiling((34 + 15 * batch_size) * 2.5)
  ))

  # ---- Thread cap for batched mode ----
  # Each CRISPRitz thread loads its own copy of the TST tree (~1.5-2 GB
  # for hg38 bMax=1). In batched mode, R's cumulative RSS grows across
  # batches despite gc(), reducing available headroom. Capping at 3
  # threads (vs 4 in non-batched mode) saves one full tree copy,
  # reducing CRISPRitz peak memory from ~6-8 GB to ~4.5-6 GB.
  # Empirical: 4 threads causes OOM at batch 9 with 14.7 GB WSL;
  # 3 threads should provide ~3 GB additional headroom.
  effective_threads <- min(as.integer(threads), 3L)
  if (!quiet && threads > 3L) {
    message("Note: Capping threads to 3 for batched indexed search ",
            "(preserves memory headroom across batches).")
  }

  # ---- Upfront CRISPRitz probe ----
  # Resolve once and reuse for all chunks. WSL can degrade after OOM kills;
  # re-probing per chunk then fails with "not found" even though the
  # installation is intact.
  crispritz_info <- locate_crispritz_bin(quiet = TRUE)
  if (is.null(crispritz_info$bin)) {
    warning("CRISPRitz executable not found. Cannot perform off-target search.")
    empty <- parse_crispritz_output("/dev/null", guide_map, pam_length = 3L,
                                    pam_side = "3prime")
    attr(empty, "searched_spacers") <- character(0)
    return(empty)
  }

  system_degraded   <- FALSE

  # ---- File-backed batch storage ----
  # Each batch's parsed results (100K-500K rows) accumulates in R memory.
  # After ~8 batches, R holds ~1 GB of hit data, reducing physical RAM
  # available to CRISPRitz (which shares memory via WSL). Spill each
  # batch to a temp RDS file and read back at the end.
  batch_tmpdir <- file.path(tempdir(), paste0("mutateR_batches_",
                                              format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(batch_tmpdir, recursive = TRUE, showWarnings = FALSE)
  batch_files <- character(0)

  searched_spacers    <- character(0)
  quarantined_spacers <- character(0)

  # ---- Pre-screen: quarantine repeat-rich guides before batch submission ----
  # In pure CRISPRitz mode, Bowtie hit counts are unavailable — use
  # linguistic complexity only. Catches Alu-derived poly-T guides (the
  # worst OOM offenders) but misses moderate-complexity repeat junctions;
  # those are caught reactively by single-guide probing.
  repeat_flags <- screen_repeat_rich(
    spacers      = unique_spacers,
    lc_threshold = 0.6,
    quiet        = quiet
  )

  if (any(repeat_flags)) {
    pre_quarantined     <- unique_spacers[repeat_flags]
    unique_spacers      <- unique_spacers[!repeat_flags]
    guide_map           <- guide_map[unique_spacers]
    quarantined_spacers <- c(quarantined_spacers, pre_quarantined)

    n_spacers <- length(unique_spacers)
    n_batches <- ceiling(n_spacers / batch_size)

    if (!quiet) message("Pre-screened ", length(pre_quarantined),
                        " repeat-rich guide(s). Remaining: ",
                        n_spacers, " guides in ", n_batches, " batches.")
  }

  # ---- Guard: all guides pre-screened out ----
  if (n_spacers == 0) {
    if (!quiet) message("All guides quarantined by pre-screening. ",
                        "No CRISPRitz search required.")
    unlink(batch_tmpdir, recursive = TRUE)
    empty <- parse_crispritz_output("/dev/null", guide_map, pam_length = 3L,
                                    pam_side = "3prime")
    attr(empty, "searched_spacers")    <- character(0)
    attr(empty, "quarantined_spacers") <- unique(quarantined_spacers)
    return(empty)
  }

  if (!quiet) {
    est_total_min <- round((34 + 15 * n_spacers) / 60, 1)
    message("Batched indexed search: ", n_spacers, " guides in ",
            n_batches, " batches of up to ", batch_size,
            " (per-batch timeout: ", batch_timeout, "s, ",
            "estimated total: ~", est_total_min, " min)")
  }

  t_start <- Sys.time()
  total_hits       <- 0L
  n_batch_ok       <- 0L
  n_batch_failed   <- 0L

  for (b in seq_len(n_batches)) {
    # ---- Degradation guard ----
    if (system_degraded) {
      if (!quiet) message("  Skipping batch ", b, "/", n_batches,
                          ": system degraded (previous probing found no working guides)")
      n_batch_failed <- n_batch_failed + 1L
      next
    }

    i_start <- (b - 1L) * batch_size + 1L
    i_end   <- min(b * batch_size, n_spacers)
    batch_spacers <- unique_spacers[i_start:i_end]
    batch_map     <- guide_map[batch_spacers]

    if (!quiet) {
      elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 1)
      message(sprintf("  Batch %d/%d: guides %d-%d (%d guides) [%s min elapsed]...",
                      b, n_batches, i_start, i_end,
                      length(batch_spacers), elapsed))
    }

    batch_error <- NULL
    batch_result <- tryCatch(
      run_crispritz_chunk(
        spacers        = batch_spacers,
        guide_map      = batch_map,
        genome_dir     = genome_dir,
        index_dir      = index_dir,
        nuclease       = nuclease,
        max_mismatches = max_mismatches,
        max_dna_bulges = max_dna_bulges,
        max_rna_bulges = max_rna_bulges,
        threads        = effective_threads,
        timeout        = batch_timeout,
        crispritz_info = crispritz_info,
        quiet          = quiet
      ),
      error = function(e) {
        if (!quiet) message("  Batch ", b, " failed: ", e$message)
        batch_error <<- e$message
        NULL
      }
    )

    # ---- Batch success ----
    if (!is.null(batch_result) && nrow(batch_result) > 0) {
      batch_rds <- file.path(batch_tmpdir, paste0("batch_", b, ".rds"))
      saveRDS(batch_result, batch_rds)
      batch_files      <- c(batch_files, batch_rds)
      n_batch_ok       <- n_batch_ok + 1L
      total_hits       <- total_hits + nrow(batch_result)
      searched_spacers <- c(searched_spacers, batch_spacers)
      rm(batch_result); gc(full = TRUE, verbose = FALSE); gc(full = TRUE, verbose = FALSE)
      Sys.sleep(3)
      reclaim_wsl_memory(quiet = quiet)
      next
    }

    # ---- OOM-specific recovery + batch retry (exit code 15 only) ----
    # OOM often reflects transient memory pressure (page cache from prior
    # batches) rather than guide-specific toxicity. Retry the full batch
    # after reclaiming memory. Two escalation levels:
    #   Level 1: drop_caches (recovers 5-7 GB page cache, ~2s)
    #   Level 2: wsl --shutdown (full VM restart, guaranteed reclamation)
    # If retries exhaust, fall through to single-guide probing which
    # identifies and quarantines the specific toxic guide(s).
    is_oom <- !is.null(batch_error) &&
      grepl("exit code 15", batch_error, fixed = TRUE)

    if (is_oom) {
      # ---- OOM: restart WSL and go straight to probing ----
      # In batched mode, OOM typically means ≥1 toxic guide in the batch.
      # Full-batch retries cannot succeed when toxic guides are present —
      # they just waste ~15 min loading the same TST tree repeatedly.
      # Restart WSL for clean memory state, then probe individually.
      if (!quiet) message("  OOM detected (exit 15). Restarting WSL for clean probing...")
      tryCatch(
        system2("wsl", args = "--shutdown", stdout = TRUE, stderr = TRUE),
        error = function(e) NULL
      )
      Sys.sleep(15)
    }


    # ---- WSL recovery after batch failure ----
    # CRISPRitz OOM (exit 15) kills the WSL process, often leaving the
    # WSL VM in a degraded state where all subsequent system2("wsl", ...)
    # calls fail immediately with exit 65535. Without recovery time,
    # single-guide probing would fail instantly for every guide, the
    # batch would be declared degraded, and remaining batches skipped.
    #
    # Pause to let the Linux OOM killer clean up, then probe WSL health.
    # If WSL hasn't recovered after two attempts, skip this batch's
    # probing entirely and let the next batch loop iteration try again.
    if (!quiet) message("  Pausing for WSL recovery after batch failure...")
    Sys.sleep(10)

    wsl_recovered <- tryCatch({
      res <- system2("wsl", args = c("bash", "-lc", shQuote("echo ok")),
                     stdout = TRUE, stderr = TRUE)
      exit <- attr(res, "status")
      if (is.null(exit)) exit <- 0L
      exit == 0L && any(grepl("ok", res))
    }, error = function(e) FALSE)

    if (!wsl_recovered) {
      if (!quiet) message("  WSL not responsive after 10s. Waiting 20s more...")
      Sys.sleep(20)
      wsl_recovered <- tryCatch({
        res <- system2("wsl", args = c("bash", "-lc", shQuote("echo ok")),
                       stdout = TRUE, stderr = TRUE)
        exit <- attr(res, "status")
        if (is.null(exit)) exit <- 0L
        exit == 0L && any(grepl("ok", res))
      }, error = function(e) FALSE)
    }

    if (!wsl_recovered) {
      if (!quiet) message("  WSL still unresponsive after 30s total. ",
                          "Skipping fallback for batch ", b,
                          ". Will retry on next batch.")
      n_batch_failed <- n_batch_failed + 1L
      next
    }

    if (!quiet) message("  WSL recovered. Proceeding with single-guide probing.")


    # ---- Batch failure: single-guide probing ----
    if (!quiet) message("  Probing batch ", b, " guides individually to isolate ",
                        "repeat-rich guides (", length(batch_spacers), " guides)...")

    probe_result <- probe_failed_batch(
      spacers        = batch_spacers,
      guide_map      = batch_map,
      genome_dir     = genome_dir,
      index_dir      = index_dir,
      nuclease       = nuclease,
      max_mismatches = max_mismatches,
      max_dna_bulges = max_dna_bulges,
      max_rna_bulges = max_rna_bulges,
      threads        = effective_threads,
      timeout        = batch_timeout,
      crispritz_info = crispritz_info,
      batch_tmpdir   = batch_tmpdir,
      batch_id       = b,
      quiet          = quiet
    )

    # Collect results
    batch_files         <- c(batch_files, probe_result$batch_files)
    searched_spacers    <- c(searched_spacers, probe_result$searched_spacers)
    quarantined_spacers <- c(quarantined_spacers, probe_result$quarantined_spacers)

    if (length(probe_result$searched_spacers) > 0) {
      n_batch_ok <- n_batch_ok + 1L
      total_hits <- total_hits + probe_result$total_hits
    } else {
      n_batch_failed <- n_batch_failed + 1L
    }

    # ---- Degradation detection ----
    # Probing quarantined ALL spacers with zero successes: system-level
    # failure rather than guide-specific.
    if (length(probe_result$quarantined_spacers) == length(batch_spacers) &&
        length(probe_result$searched_spacers) == 0) {
      system_degraded <- TRUE
      if (!quiet) message("  System degradation detected: all ",
                          length(batch_spacers), " probed guides failed. ",
                          "Skipping remaining batches.")
    }
  }

  # ---- Concatenate all batch results from disk ----
  if (length(batch_files) == 0) {
    t_elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 1)
    warning("All CRISPRitz batches failed (", n_batches, " batches, ",
            t_elapsed, " min). No off-target hits recovered.")

    unlink(batch_tmpdir, recursive = TRUE)
    empty <- parse_crispritz_output("/dev/null", guide_map, pam_length = 3L,
                                    pam_side = "3prime")
    attr(empty, "searched_spacers") <- unique(searched_spacers)
    return(empty)
  }

  if (!quiet) message("Reassembling results from ", length(batch_files),
                      " batch files...")
  batch_dfs <- lapply(batch_files, readRDS)
  combined  <- do.call(rbind, batch_dfs)
  rm(batch_dfs); gc(verbose = FALSE)
  unlink(batch_tmpdir, recursive = TRUE)

  t_elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 1)
  if (!quiet) {
    message(sprintf(
      "Batched search complete: %s total hits from %d/%d batches in %s min",
      format(total_hits, big.mark = ","), n_batch_ok, n_batches, t_elapsed
    ))
    if (n_batch_failed > 0) {
      message("  WARNING: ", n_batch_failed, " batch(es) failed entirely.")
    }
    if (length(quarantined_spacers) > 0) {
      message("  Quarantined (repeat-rich): ", length(quarantined_spacers),
              " guide(s) assigned specificity = 0.")
    }
  }

  attr(combined, "searched_spacers")    <- unique(searched_spacers)
  attr(combined, "quarantined_spacers") <- unique(quarantined_spacers)
  return(combined)
}


#' Auto-detect on-target score cutoff from scoring method metadata
#'
#' Replicates the cutoff logic from \code{\link{assemble_grna_pairs}} so that
#' \code{\link{refine_with_bulges}} can re-evaluate the \code{recommended} flag
#' without requiring the cutoff to be threaded from the caller.
#'
#' @param scored_grnas GRanges with \code{scoring_method} metadata column.
#'
#' @return Numeric. The appropriate on-target score cutoff.
#' @keywords internal
detect_score_cutoff <- function(scored_grnas) {
  method <- if ("scoring_method" %in% names(mcols(scored_grnas))) {
    unique(mcols(scored_grnas)$scoring_method)[1]
  } else {
    "ruleset1"
  }

  regression_models <- c("deepcpf1", "deepspcas9")
  zscore_models     <- c("ruleset3")

  if (!is.na(method) && tolower(method) %in% regression_models) return(50)
  if (!is.na(method) && tolower(method) %in% zscore_models)     return(0.1)
  return(0.5)
}


#' Refine specificity scores for recommended-pair gRNAs via CRISPRitz bulge search
#'
#' Called by \code{\link{run_mutateR}} in hybrid mode after pair assembly.
#' Extracts unique protospacers from recommended pairs, runs a comprehensive
#' CRISPRitz search (mismatches + bulges) on that small subset, and updates
#' specificity scores in both the pairs table and the scored gRNA GRanges.
#'
#' The comprehensive CRISPRitz results \strong{replace} (not supplement) the
#' initial Bowtie mismatch-only scores for the refined gRNAs, since CRISPRitz
#' in bulge mode returns all hit types including mismatch-only hits. This
#' avoids double-counting.
#'
#' @param pairs_df data.frame. From \code{\link{assemble_grna_pairs}}.
#'        Must contain columns: \code{protospacer_sequence_5p}, \code{cut_site_5p},
#'        \code{protospacer_sequence_3p}, \code{cut_site_3p}, \code{recommended},
#'        \code{ontarget_score_5p}, \code{ontarget_score_3p},
#'        \code{specificity_score_5p}, \code{specificity_score_3p}.
#' @param scored_grnas GRanges. Full scored gRNA object (all guides, not just
#'        recommended). Used for gRNA metadata lookup and on-target exclusion.
#' @param genome BSgenome object.
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a".
#' @param exon_gr GRanges. Exon coordinates for off-target annotation.
#' @param max_mismatches Integer. Maximum mismatches for CRISPRitz search.
#' @param max_dna_bulges Integer. Maximum DNA bulge size.
#' @param max_rna_bulges Integer. Maximum RNA bulge size.
#' @param bulge_penalty Numeric. Per-bulge-nucleotide CFD penalty.
#' @param threads Integer. CRISPRitz thread count.
#' @param timeout Integer. CRISPRitz timeout in seconds.
#' @param min_pair_specificity Numeric. Minimum pair-level specificity (default
#'   10) for recommended status, using the additive off-target burden model.
#'   See \code{\link{assemble_grna_pairs}} and \code{\link{compute_pair_specificity}}.
#' @param score_cutoff Numeric or NULL. On-target score cutoff for recommended
#'        status. If NULL (default), auto-detected from \code{scored_grnas} metadata.
#' @param cache_dir Character or NULL. Cache directory.
#' @param detail_level Character. Controls retention of bulge-pass off-target
#'   hit details (default \code{"compact"}). Passed through from
#'   \code{\link{score_offtargets}} / \code{\link{run_mutateR}}. See
#'   \code{\link{score_offtargets}} for level descriptions. Bulge-pass
#'   details are attached to \code{scored_grnas} as
#'   \code{attr(, "bulge_offtarget_details")}.
#' @param quiet Logical.
#'
#' @return A named list:
#'   \describe{
#'     \item{pairs}{Updated pairs_df with refined specificity scores and
#'           re-evaluated \code{recommended} flag.}
#'     \item{scored_grnas}{Updated GRanges with refined specificity scores
#'           for the subset of gRNAs that were re-searched.}
#'   }
#' @keywords internal
refine_with_bulges <- function(pairs_df,
                               scored_grnas,
                               genome,
                               nuclease,
                               exon_gr,
                               max_mismatches,
                               max_dna_bulges,
                               max_rna_bulges,
                               bulge_penalty = 0.2,
                               threads = 4L,
                               timeout = 600L,
                               min_pair_specificity = 10,
                               score_cutoff = NULL,
                               cache_dir = NULL,
                               detail_level = c("compact", "full", "file", "none"),
                               quiet = FALSE) {

  # ---- 0. Early-exit guard ----
  unchanged <- list(pairs = pairs_df, scored_grnas = scored_grnas)

  if (is.null(pairs_df) || !is.data.frame(pairs_df) || nrow(pairs_df) == 0) {
    return(unchanged)
  }

  rec <- pairs_df[pairs_df$recommended == TRUE, , drop = FALSE]
  if (nrow(rec) == 0) {
    if (!quiet) message("No recommended pairs to refine. Skipping bulge refinement.")
    return(unchanged)
  }

  detail_level <- match.arg(detail_level)

  # ---- 1. Check CRISPRitz availability ----
  if (!check_mutater_env()) {
    warning("mutateR Python environment not available for CRISPRitz bulge refinement.\n",
            "Skipping bulge pass. Specificity scores reflect mismatch-only (Bowtie) analysis.\n",
            "To enable bulge refinement, run install_mutater_env() and restart R.")
    return(unchanged)
  }

  # ---- 2. Extract target protospacers from recommended pairs ----
  target_spacers <- unique(toupper(c(
    as.character(rec$protospacer_sequence_5p),
    as.character(rec$protospacer_sequence_3p)
  )))

  if (!quiet) message("\n--- Hybrid mode: bulge refinement ---")
  if (!quiet) message("Target protospacers: ", length(target_spacers),
                      " (from ", nrow(rec), " recommended pairs)")

  # ---- 3. Find matching gRNAs in scored_grnas ----
  # Match on protospacer sequence to find all gRNAs sharing these spacers
  # (a non-recommended pair may share a protospacer with a recommended one)
  all_spacers <- toupper(as.character(mcols(scored_grnas)$protospacer_sequence))
  subset_idx  <- which(all_spacers %in% target_spacers)

  if (length(subset_idx) == 0) {
    warning("No matching gRNAs found in scored_grnas for bulge refinement. ",
            "This should not happen — check protospacer_sequence consistency.")
    return(unchanged)
  }

  subset_gr <- scored_grnas[subset_idx]

  if (!quiet) message("Matched ", length(subset_gr), " gRNAs in scored_grnas ",
                      "(", length(target_spacers), " unique protospacers)")

  # ---- 3b. Pre-screen repeat-rich guides before CRISPRitz submission ----
  # In hybrid mode, Phase 1 Bowtie scores are available as n_offtargets.
  # Guides exceeding the Bowtie threshold or with low linguistic complexity
  # are removed from the CRISPRitz refinement set. They retain their
  # Phase 1 Bowtie specificity scores — conservative (mismatch-only
  # underestimates off-target burden) but still informative, and far
  # better than crashing CRISPRitz.
  subset_spacers   <- toupper(as.character(mcols(subset_gr)$protospacer_sequence))
  subset_ot_counts <- if ("n_offtargets" %in% names(mcols(subset_gr))) {
    mcols(subset_gr)$n_offtargets
  } else {
    NULL
  }

  repeat_flags <- screen_repeat_rich(
    spacers           = subset_spacers,
    bowtie_hit_counts = subset_ot_counts,
    bowtie_threshold  = 5000L,
    lc_threshold      = 0.6,
    quiet             = quiet
  )

  if (any(repeat_flags)) {
    n_pre_screened <- sum(repeat_flags)
    if (!quiet) message("Pre-screened ", n_pre_screened,
                        " repeat-rich guide(s) from bulge refinement set. ",
                        "These retain Bowtie-only specificity scores.")

    subset_idx <- subset_idx[!repeat_flags]
    subset_gr  <- scored_grnas[subset_idx]
    target_spacers <- unique(toupper(as.character(
      mcols(subset_gr)$protospacer_sequence
    )))
  }

  # Re-check: any guides left to refine?
  if (length(subset_idx) == 0) {
    if (!quiet) message("No guides remaining for bulge refinement after pre-screening.")
    return(unchanged)
  }

  # ---- 4. Deduplicate protospacers for search ----
  dedup <- deduplicate_protospacers(subset_gr)

  if (!quiet) {
    est_seconds <- length(dedup$unique_spacers) * max_mismatches * 0.5
    est_minutes <- round(est_seconds / 60, 1)
    message("Launching CRISPRitz search: ", length(dedup$unique_spacers),
            " unique spacers, mm=", max_mismatches,
            ", bDNA=", max_dna_bulges, ", bRNA=", max_rna_bulges,
            " (est. ~", est_minutes, " min)")
  }

  # ---- 5. Setup cache ----
  if (is.null(cache_dir)) {
    cache_dir <- tools::R_user_dir("mutateR", "cache")
  }

  # ---- 6a. Run CRISPRitz backend on subset ----
  t_start <- Sys.time()

  hit_df <- tryCatch(
    run_crispritz_backend(
      unique_spacers = dedup$unique_spacers,
      guide_map      = dedup$guide_map,
      genome         = genome,
      nuclease       = nuclease,
      max_mismatches = as.integer(max_mismatches),
      max_dna_bulges = as.integer(max_dna_bulges),
      max_rna_bulges = as.integer(max_rna_bulges),
      threads        = as.integer(threads),
      timeout        = as.integer(timeout),
      cache_dir      = cache_dir,
      quiet          = quiet
    ),
    error = function(e) {
      if (!quiet) message("Batch bulge search failed: ", e$message)
      NULL
    }
  )

  # ---- 6b. Single-guide sequential fallback ----
  # CRISPRitz indexed bulge search can crash for multi-guide batches
  # (known issue at bMax>=2). Fall back to processing guides individually.
  if (is.null(hit_df) && length(dedup$unique_spacers) > 1) {
    if (!quiet) message("Attempting single-guide sequential fallback ",
                        "(", length(dedup$unique_spacers), " guides)...")

    hit_dfs <- vector("list", length(dedup$unique_spacers))
    n_succeeded <- 0L
    n_failed    <- 0L

    for (i in seq_along(dedup$unique_spacers)) {
      spacer_i <- dedup$unique_spacers[i]
      map_i    <- dedup$guide_map[spacer_i]

      single_hit <- tryCatch(
        run_crispritz_backend(
          unique_spacers = spacer_i,
          guide_map      = map_i,
          genome         = genome,
          nuclease       = nuclease,
          max_mismatches = as.integer(max_mismatches),
          max_dna_bulges = as.integer(max_dna_bulges),
          max_rna_bulges = as.integer(max_rna_bulges),
          threads        = min(as.integer(threads), 4L),
          timeout        = as.integer(timeout),
          cache_dir      = cache_dir,
          quiet          = TRUE
        ),
        error = function(e) NULL
      )

      if (!is.null(single_hit) && nrow(single_hit) > 0) {
        hit_dfs[[i]] <- single_hit
        n_succeeded  <- n_succeeded + 1L
      } else {
        n_failed <- n_failed + 1L
      }

      if (!quiet && (i %% 5 == 0 || i == length(dedup$unique_spacers))) {
        message("  Progress: ", i, "/", length(dedup$unique_spacers),
                " (", n_succeeded, " succeeded, ", n_failed, " failed)")
      }
    }

    non_null <- Filter(function(x) !is.null(x), hit_dfs)
    if (length(non_null) > 0) {
      hit_df <- do.call(rbind, non_null)
      if (!quiet) message("Sequential fallback complete: ",
                          n_succeeded, "/", length(dedup$unique_spacers),
                          " guides succeeded, ", nrow(hit_df), " total hits.")
    } else {
      if (!quiet) message("All single-guide searches failed.")
    }
  }

  if (is.null(hit_df)) {
    warning("CRISPRitz bulge search failed (batch and sequential). ",
            "Retaining mismatch-only specificity scores.")
    return(unchanged)
  }

  t_elapsed <- difftime(Sys.time(), t_start, units = "secs")
  if (!quiet) message("CRISPRitz bulge search completed in ",
                      round(as.numeric(t_elapsed), 1), " seconds. ",
                      nrow(hit_df), " hits found.")

  # ---- 7a. Score and aggregate on subset ----
  agg <- score_and_aggregate(
    hit_df        = hit_df,
    grna_gr       = subset_gr,
    exon_gr       = exon_gr,
    nuclease      = nuclease,
    bulge_penalty = bulge_penalty
  )

  # ---- 7b. Handle bulge hit details per detail_level ----
  bulge_details      <- NULL
  bulge_details_path <- NULL

  if (detail_level == "full") {
    bulge_details <- agg$hit_details

  } else if (detail_level == "compact") {
    bulge_details <- compact_hit_details(agg$hit_details)

  } else if (detail_level == "file") {
    if (is.null(cache_dir)) cache_dir <- tools::R_user_dir("mutateR", "cache")
    bulge_details_path <- file.path(
      cache_dir,
      paste0("bulge_offtarget_details_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    )
    saveRDS(agg$hit_details, bulge_details_path, compress = "gzip")
    if (!quiet) message("Bulge off-target details written to: ", bulge_details_path)
    bulge_details <- compact_hit_details(agg$hit_details)
  }
  # "none" → bulge_details stays NULL

  attr(scored_grnas, "bulge_offtarget_details")      <- bulge_details
  attr(scored_grnas, "bulge_offtarget_details_path") <- bulge_details_path

  # ---- 7c. Defensive: ensure OT columns exist in scored_grnas ----
  # Phase 1 may have failed (env unavailable, timeout, transient crash),
  # leaving scored_grnas without specificity columns. Initialise them
  # so the indexed assignment in Step 9 doesn't crash with a length
  # mismatch when writing into a NULL column.
  if (!"specificity_score" %in% names(mcols(scored_grnas))) {
    mcols(scored_grnas)$specificity_score <- NA_real_
  }
  if (!"n_offtargets" %in% names(mcols(scored_grnas))) {
    mcols(scored_grnas)$n_offtargets <- NA_integer_
  }
  if (!"has_exonic_offtarget" %in% names(mcols(scored_grnas))) {
    mcols(scored_grnas)$has_exonic_offtarget <- FALSE
  }

  # ---- 8. Build lookup keyed by (protospacer, cut_site) ----
  # This handles the case where the same protospacer appears at different
  # genomic loci (different exons), which have different on-target exclusion
  # windows and therefore different specificity scores.
  lookup <- data.frame(
    protospacer          = toupper(as.character(mcols(subset_gr)$protospacer_sequence)),
    cut_site             = mcols(subset_gr)$cut_site,
    specificity_score    = agg$per_grna$specificity_score,
    n_offtargets         = agg$per_grna$n_offtargets,
    has_exonic_offtarget = agg$per_grna$has_exonic_offtarget,
    stringsAsFactors     = FALSE
  )
  lookup$key <- paste(lookup$protospacer, lookup$cut_site, sep = ":")

  # ---- 9. Update scored_grnas for refined subset ----
  mcols(scored_grnas)$specificity_score[subset_idx]    <- agg$per_grna$specificity_score
  mcols(scored_grnas)$n_offtargets[subset_idx]         <- agg$per_grna$n_offtargets
  mcols(scored_grnas)$has_exonic_offtarget[subset_idx] <- agg$per_grna$has_exonic_offtarget

  # ---- 10. Update pairs_df via protospacer + cut_site matching ----
  # 5' gRNA
  key_5p   <- paste(toupper(pairs_df$protospacer_sequence_5p),
                    pairs_df$cut_site_5p, sep = ":")
  match_5p <- match(key_5p, lookup$key)
  has_5p   <- !is.na(match_5p)

  if (any(has_5p)) {
    pairs_df$specificity_score_5p[has_5p] <- lookup$specificity_score[match_5p[has_5p]]
    pairs_df$n_offtargets_5p[has_5p]      <- lookup$n_offtargets[match_5p[has_5p]]
  }

  # 3' gRNA
  key_3p   <- paste(toupper(pairs_df$protospacer_sequence_3p),
                    pairs_df$cut_site_3p, sep = ":")
  match_3p <- match(key_3p, lookup$key)
  has_3p   <- !is.na(match_3p)

  if (any(has_3p)) {
    pairs_df$specificity_score_3p[has_3p] <- lookup$specificity_score[match_3p[has_3p]]
    pairs_df$n_offtargets_3p[has_3p]      <- lookup$n_offtargets[match_3p[has_3p]]
  }

  # ---- 11. Recalculate pair_specificity and re-evaluate recommended ----
  pairs_df$pair_specificity <- compute_pair_specificity(
    pairs_df$specificity_score_5p,
    pairs_df$specificity_score_3p
  )

  # Determine score_cutoff if not provided
  if (is.null(score_cutoff)) {
    score_cutoff <- detect_score_cutoff(scored_grnas)
  }

  n_rec_before <- sum(pairs_df$recommended, na.rm = TRUE)

  pairs_df$recommended <- with(pairs_df,
                               !is.na(ontarget_score_5p) & ontarget_score_5p >= score_cutoff &
                                 !is.na(ontarget_score_3p) & ontarget_score_3p >= score_cutoff &
                                 (is.na(pair_specificity) | pair_specificity >= min_pair_specificity))

  n_rec_after <- sum(pairs_df$recommended, na.rm = TRUE)

  # ---- 12. Report ----
  if (!quiet) {
    n_lost    <- n_rec_before - n_rec_after
    n_refined <- sum(has_5p | has_3p)

    message("\n--- Bulge refinement summary ---")
    message("  Pairs with updated scores: ", n_refined, " / ", nrow(pairs_df))
    message("  Recommended before bulge:  ", n_rec_before)
    message("  Recommended after bulge:   ", n_rec_after)
    if (n_lost > 0) {
      message("  WARNING: ", n_lost, " pair(s) lost recommended status ",
              "due to bulge off-targets reducing pair specificity below ", min_pair_specificity)
    }
    if (n_rec_after == 0) {
      message("  WARNING: No recommended pairs remain after bulge refinement. ",
              "Consider relaxing min_pair_specificity or reviewing off-target landscape.")
    }
    message("--------------------------------\n")
  }

  return(list(pairs = pairs_df, scored_grnas = scored_grnas))
}

#' Score off-target hits (CFD), annotate genomic context, and aggregate per gRNA
#'
#' Backend-agnostic function that takes a standardised hit data.frame (from
#' any search backend), computes per-hit CFD scores, annotates exonic overlaps,
#' and aggregates to per-gRNA MIT-style specificity scores.
#'
#' Performance: all aggregation is vectorised via index expansion + tapply.
#' Scales linearly with total (hit x gRNA) pairs.
#'
#' @param hit_df data.frame. Standardised hit table with columns: guide_seq,
#'        grna_indices, chr, pos, strand, offtarget_seq, pam_gen, n_mismatches,
#'        bulge_type, bulge_size.
#' @param grna_gr GRanges. Original gRNA object (for on-target exclusion).
#' @param exon_gr GRanges. Exon coordinates (for exonic annotation).
#' @param nuclease Character. One of "Cas9", "Cas12a", or "enCas12a".
#' @param bulge_penalty Numeric. Penalty per bulge nucleotide.
#'
#' @return A named list:
#'   \describe{
#'     \item{per_grna}{data.frame with columns: grna_index, specificity_score,
#'           n_offtargets, has_exonic_offtarget.}
#'     \item{hit_details}{data.frame of all scored and annotated off-target hits.}
#'   }
#' @keywords internal
score_and_aggregate <- function(hit_df, grna_gr, exon_gr, nuclease, bulge_penalty = 0.2) {

  # ---- 0. Handle empty input ----
  n_grnas <- length(grna_gr)

  empty_per_grna <- data.frame(
    grna_index           = seq_len(n_grnas),
    specificity_score    = rep(100, n_grnas),
    n_offtargets         = rep(0L, n_grnas),
    has_exonic_offtarget = rep(FALSE, n_grnas),
    stringsAsFactors     = FALSE
  )

  if (nrow(hit_df) == 0) {
    return(list(
      per_grna    = empty_per_grna,
      hit_details = hit_df
    ))
  }

  # ---- 1. Per-hit CFD scoring ----
  guide_seqs <- hit_df$guide_seq
  ot_seqs    <- hit_df$offtarget_seq

  hit_df$cfd_score <- NA_real_

  is_mismatch_only <- hit_df$bulge_type == "X"
  is_bulge         <- !is_mismatch_only

  # ---- 1a. Mismatch-only hits: batch CFD scoring ----
  if (any(is_mismatch_only)) {

    mm_guides  <- guide_seqs[is_mismatch_only]
    mm_targets <- ot_seqs[is_mismatch_only]

    mm_cfd <- tryCatch({
      valid_len <- nchar(mm_guides) == nchar(mm_targets)
      scores <- rep(NA_real_, length(mm_guides))

      if (any(valid_len)) {
        guide_length <- nchar(mm_guides[valid_len][1])

        if (guide_length == 20 && nuclease == "Cas9") {
          mm_pams <- hit_df$pam_gen[is_mismatch_only][valid_len]
          if (all(is.na(mm_pams))) mm_pams <- rep("NGG", sum(valid_len))
          mm_pams[is.na(mm_pams)] <- "NGG"

          cfd_result <- crisprScore::getCFDScores(
            spacers      = mm_guides[valid_len],
            protospacers = mm_targets[valid_len],
            pams         = mm_pams
          )
          scores[valid_len] <- if (is.data.frame(cfd_result)) {
            as.numeric(cfd_result$score)
          } else {
            as.numeric(cfd_result)
          }
        } else {
          scores[valid_len] <- 0.5 ^ hit_df$n_mismatches[is_mismatch_only][valid_len]
        }
      }

      scores
    }, error = function(e) {
      warning("CFD scoring failed for mismatch-only hits: ", e$message,
              "\nFalling back to simple mismatch-count heuristic.")
      0.5 ^ hit_df$n_mismatches[is_mismatch_only]
    })

    hit_df$cfd_score[is_mismatch_only] <- mm_cfd
  }

  # ---- 1b. Bulge-containing hits: batch where possible ----
  if (any(is_bulge)) {

    bulge_idx     <- which(is_bulge)
    bulge_guides  <- guide_seqs[bulge_idx]
    bulge_targets <- ot_seqs[bulge_idx]
    bulge_sizes   <- hit_df$bulge_size[bulge_idx]
    bulge_mm      <- hit_df$n_mismatches[bulge_idx]

    bulge_cfd <- tryCatch({
      g_len <- nchar(bulge_guides)
      t_len <- nchar(bulge_targets)
      compare_len <- pmin(g_len, t_len)

      g_trim <- substr(bulge_guides, 1, compare_len)
      t_trim <- substr(bulge_targets, 1, compare_len)

      is_cas9_20 <- (compare_len == 20) & (nuclease == "Cas9")

      mismatch_cfd <- rep(NA_real_, length(bulge_idx))

      if (any(is_cas9_20)) {
        batch_result <- tryCatch({
          res <- crisprScore::getCFDScores(
            spacers      = g_trim[is_cas9_20],
            protospacers = t_trim[is_cas9_20]
          )
          if (is.data.frame(res)) as.numeric(res$score) else as.numeric(res)
        }, error = function(e) {
          0.5 ^ bulge_mm[is_cas9_20]
        })
        mismatch_cfd[is_cas9_20] <- batch_result
      }

      if (any(!is_cas9_20)) {
        mismatch_cfd[!is_cas9_20] <- 0.5 ^ bulge_mm[!is_cas9_20]
      }

      mismatch_cfd * (bulge_penalty ^ bulge_sizes)

    }, error = function(e) {
      warning("CFD scoring failed for bulge hits: ", e$message,
              "\nUsing bulge penalty only.")
      bulge_penalty ^ bulge_sizes
    })

    hit_df$cfd_score[bulge_idx] <- bulge_cfd
  }

  hit_df$cfd_score <- pmin(pmax(hit_df$cfd_score, 0), 1)
  hit_df$cfd_score[is.na(hit_df$cfd_score)] <- 0

  # ---- 2. Genomic annotation (vectorised) ----
  hit_gr <- tryCatch({
    adjusted_width <- nchar(hit_df$offtarget_seq) +
      ifelse(hit_df$bulge_type == "RNA", hit_df$bulge_size, 0L)

    gr <- GenomicRanges::GRanges(
      seqnames = hit_df$chr,
      ranges   = IRanges::IRanges(start = hit_df$pos,
                                  width = adjusted_width),
      strand   = hit_df$strand
    )

    tryCatch(
      GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(exon_gr)[1],
      error = function(e) NULL
    )

    gr
  }, error = function(e) {
    warning("Could not build GRanges for off-target annotation: ", e$message)
    NULL
  })

  hit_df$is_exonic <- FALSE
  hit_df$exon_hit  <- NA_character_

  if (!is.null(hit_gr) && length(hit_gr) == nrow(hit_df)) {
    overlaps <- tryCatch(
      GenomicRanges::findOverlaps(hit_gr, exon_gr, type = "any"),
      error = function(e) {
        warning("Exon overlap detection failed: ", e$message)
        NULL
      }
    )

    if (!is.null(overlaps) && length(overlaps) > 0) {
      q_hits <- S4Vectors::queryHits(overlaps)
      s_hits <- S4Vectors::subjectHits(overlaps)

      exonic_hit_idx <- unique(q_hits)
      hit_df$is_exonic[exonic_hit_idx] <- TRUE

      exon_ranks <- mcols(exon_gr)$rank[s_hits]
      exon_hit_strings <- tapply(exon_ranks, q_hits, function(x) {
        paste(unique(x), collapse = ",")
      })
      hit_df$exon_hit[as.integer(names(exon_hit_strings))] <- as.character(exon_hit_strings)
    }
  }

  # ---- 3. Per-gRNA aggregation (vectorised) ----
  grna_chr <- as.character(GenomicRanges::seqnames(grna_gr))
  grna_cut <- mcols(grna_gr)$cut_site

  grna_chr_style <- tryCatch(
    GenomeInfoDb::seqlevelsStyle(grna_gr)[1],
    error = function(e) "UCSC"
  )

  hit_chr_mapped <- tryCatch({
    unique_chrs <- unique(hit_df$chr)
    if (length(unique_chrs) == 0) return(character(0))
    tmp_gr <- GenomicRanges::GRanges(
      seqnames = unique_chrs,
      ranges   = IRanges::IRanges(start = 1, width = 1)
    )
    GenomeInfoDb::seqlevelsStyle(tmp_gr) <- grna_chr_style
    mapping <- setNames(as.character(GenomicRanges::seqnames(tmp_gr)),
                        unique_chrs)
    mapping[hit_df$chr]
  }, error = function(e) {
    hit_df$chr
  })

  # 3a. Expand (hit, gRNA) pairs
  idx_lengths  <- lengths(hit_df$grna_indices)
  exp_hit_idx  <- rep(seq_len(nrow(hit_df)), idx_lengths)
  exp_grna_idx <- unlist(hit_df$grna_indices, use.names = FALSE)

  if (length(exp_grna_idx) == 0) {
    return(list(
      per_grna    = empty_per_grna,
      hit_details = hit_df
    ))
  }

  # 3b. Pull hit-level and gRNA-level metadata via index
  exp_chr    <- hit_chr_mapped[exp_hit_idx]
  exp_pos    <- hit_df$pos[exp_hit_idx]
  exp_cfd    <- hit_df$cfd_score[exp_hit_idx]
  exp_exonic <- hit_df$is_exonic[exp_hit_idx]

  exp_grna_chr <- grna_chr[exp_grna_idx]
  exp_grna_cut <- grna_cut[exp_grna_idx]

  # 3c. On-target exclusion
  is_ontarget <- (exp_chr == exp_grna_chr) & (abs(exp_pos - exp_grna_cut) <= 10)
  keep <- !is_ontarget

  # 3d. Subset to off-target pairs
  ot_grna_idx <- exp_grna_idx[keep]
  ot_cfd      <- exp_cfd[keep]
  ot_exonic   <- exp_exonic[keep]

  # 3e. Aggregate per gRNA
  if (length(ot_grna_idx) > 0) {
    gi_factor <- factor(ot_grna_idx, levels = seq_len(n_grnas))

    cfd_sums    <- as.numeric(tapply(ot_cfd, gi_factor, sum))
    n_ot_counts <- as.integer(tapply(ot_cfd > 0.01, gi_factor, sum))
    has_exonic  <- as.logical(tapply(ot_exonic & ot_cfd > 0.1, gi_factor, any))

    cfd_sums[is.na(cfd_sums)]       <- 0
    n_ot_counts[is.na(n_ot_counts)] <- 0L
    has_exonic[is.na(has_exonic)]    <- FALSE
  } else {
    cfd_sums    <- rep(0, n_grnas)
    n_ot_counts <- rep(0L, n_grnas)
    has_exonic  <- rep(FALSE, n_grnas)
  }

  spec_scores <- 100 / (1 + cfd_sums)

  # ---- 4. Assemble results ----
  per_grna <- data.frame(
    grna_index           = seq_len(n_grnas),
    specificity_score    = round(spec_scores, 2),
    n_offtargets         = n_ot_counts,
    has_exonic_offtarget = has_exonic,
    stringsAsFactors     = FALSE
  )

  # ---- 5. Diagnostic: mismatch vs bulge CFD breakdown ----
  # Informational only — no scoring behaviour change. Helps diagnose
  # specificity deflation when bulge hits dominate the CFD sum
  # (observed in INS all-CRISPRitz: specificity 0.45-7.62).
  if (nrow(hit_df) > 0) {
    is_mm    <- hit_df$bulge_type == "X"
    n_mm     <- sum(is_mm)
    n_bulge  <- sum(!is_mm)

    mm_cfd_sum    <- sum(hit_df$cfd_score[is_mm], na.rm = TRUE)
    bulge_cfd_sum <- sum(hit_df$cfd_score[!is_mm], na.rm = TRUE)
    total_cfd_sum <- mm_cfd_sum + bulge_cfd_sum

    message("  Off-target hit breakdown:")
    message("    Mismatch-only: ", format(n_mm, big.mark = ","),
            " hits, CFD sum = ", round(mm_cfd_sum, 2),
            " (", round(100 * mm_cfd_sum / max(total_cfd_sum, 1e-9), 1), "%)")
    message("    Bulge hits:    ", format(n_bulge, big.mark = ","),
            " hits, CFD sum = ", round(bulge_cfd_sum, 2),
            " (", round(100 * bulge_cfd_sum / max(total_cfd_sum, 1e-9), 1), "%)")
    message("    Bulge penalty: ", bulge_penalty, " per nucleotide")

    if (n_bulge > 0) {
      mean_bulge_cfd <- round(bulge_cfd_sum / n_bulge, 4)
      message("    Mean per-bulge-hit CFD: ", mean_bulge_cfd)
    }

    if (total_cfd_sum > 0 && bulge_cfd_sum > mm_cfd_sum * 2) {
      message("    NOTE: Bulge hits dominate CFD sum (>2x mismatch contribution). ",
              "Consider comparing Bowtie-only vs CRISPRitz scores for same guides ",
              "to assess whether bulge_penalty (", bulge_penalty, ") needs adjustment.")
    }

    # Per-gRNA summary: median and range of specificity
    valid_spec <- per_grna$specificity_score[per_grna$specificity_score < 100]
    if (length(valid_spec) > 0) {
      message("    Specificity (guides with OT hits): median = ",
              round(median(valid_spec), 1),
              ", range = [", round(min(valid_spec), 2),
              ", ", round(max(valid_spec), 2), "]")
    }
  }

  return(list(
    per_grna    = per_grna,
    hit_details = hit_df
  ))
}

#' Compact off-target hit details for memory-efficient storage
#'
#' Reduces a full off-target hit table to a manageable size by filtering out
#' negligible-CFD hits and capping the number of retained hits per gRNA.
#' Called by \code{\link{score_offtargets}} when \code{detail_level} is
#' \code{"compact"} or \code{"file"}.
#'
#' Aggregation accuracy is unaffected: per-gRNA specificity scores are always
#' computed from the \emph{complete} hit set before compaction.
#'
#' @param hit_details data.frame. Full scored hit table as returned by
#'   \code{score_and_aggregate()$hit_details}. Must contain columns
#'   \code{cfd_score} and list-column \code{grna_indices}.
#' @param max_per_grna Integer. Maximum hits retained per gRNA, ranked by
#'   descending CFD score (default 200).
#' @param min_cfd Numeric. Minimum CFD score for a hit to be retained
#'   (default 0.01). Hits below this threshold contribute negligibly to
#'   specificity and are discarded before ranking.
#'
#' @return A filtered data.frame with the same column schema as the input.
#'   A hit row is retained if it ranks in the top \code{max_per_grna} for
#'   \emph{any} of its associated gRNAs. The \code{grna_indices} list-column
#'   is preserved unmodified.
#'
#' @keywords internal
compact_hit_details <- function(hit_details,
                                max_per_grna = 200L,
                                min_cfd = 0.01) {

  if (is.null(hit_details) || nrow(hit_details) == 0) {
    return(hit_details)
  }

  # ---- 1. Filter to hits above CFD threshold ----
  above_thresh <- which(hit_details$cfd_score > min_cfd)
  if (length(above_thresh) == 0) {
    return(hit_details[integer(0), , drop = FALSE])
  }
  filtered <- hit_details[above_thresh, , drop = FALSE]

  # ---- 2. Expand (hit_row, grna_index) pairs ----
  idx_lengths  <- lengths(filtered$grna_indices)
  exp_hit_row  <- rep(seq_len(nrow(filtered)), idx_lengths)
  exp_grna_idx <- unlist(filtered$grna_indices, use.names = FALSE)
  exp_cfd      <- filtered$cfd_score[exp_hit_row]

  if (length(exp_grna_idx) == 0) {
    return(filtered)
  }

  # ---- 3. Rank by CFD within each gRNA, keep top N ----
  ord <- order(exp_grna_idx, -exp_cfd)
  exp_grna_sorted <- exp_grna_idx[ord]
  exp_hit_sorted  <- exp_hit_row[ord]

  # Within-group rank via run-length encoding on sorted gRNA index
  grna_runs <- rle(exp_grna_sorted)
  ranks <- unlist(lapply(grna_runs$lengths, seq_len), use.names = FALSE)

  keep_expanded <- which(ranks <= max_per_grna)

  # ---- 4. Recover unique hit rows that survived ----
  keep_rows <- sort(unique(exp_hit_sorted[keep_expanded]))

  # ---- 5. Return filtered subset ----
  filtered[keep_rows, , drop = FALSE]
}
