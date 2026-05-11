#' @title Run the complete mutateR workflow
#'
#' @description Executes all major steps of the mutateR pipeline end‑to‑end:
#' 1. Gene/Transcript retrieval
#' 2. Exon phase mapping
#' 3. gRNA finding & scoring (Cas9, Cas12a, or enCas12a)
#' 4. Exon-deletion pair assembly (Auto-detects scoring thresholds unless overridden)
#' 5. Genotyping primer design (Runs only on recommended pairs)
#' 6. Visualisation
#'
#' @param gene_id Character. Gene symbol or Ensembl Gene ID (ENSG...).
#' @param species Character. e.g. "hsapiens", "mmusculus", "drerio".
#' @param genome  BSgenome object (e.g. BSgenome.Hsapiens.UCSC.hg38).
#' @param nuclease Either a character string ("Cas9" (NGG PAM), "Cas12a"
#'        (TTTV PAM), or "enCas12a" (TTTN PAM); defaults to "Cas9") OR a
#'        \code{NucleaseSpec} object built with \code{\link{nuclease_spec}}.
#'        When a custom (non-canonical) \code{NucleaseSpec} is supplied,
#'        on-target and off-target scoring are automatically skipped: the
#'        pipeline runs site-finding, phase filtering, pair assembly, primer
#'        design, and plotting only. For lightweight site-scanning with a
#'        custom nuclease, see also \code{\link{run_mutateR_custom}}.
#' @param transcript_id Optional Ensembl transcript ID to override canonical.
#' @param score_method Character. On‑target scoring model. If NULL (default), autoselects:
#'        - Cas9: "ruleset1" (alternatives: "deephf", "deepspcas9", "ruleset3")
#'        - Cas12a: "deepcpf1"
#'        - enCas12a: "enpamgb"
#' @param tracr Character. For Rule Set 3 scoring - one of "Chen2013" (default) or "Hsu2013".
#' @param deephf_var Character. For DeepHF scoring - one of "wt", "wt_u6" (default), "wt_t7", "esp", or "hf".
#'                See \code{\link{recommend_deephf_model}} for guidance on model selection based on
#'                experimental design (delivery method, Cas9 variant).
#' @param min_score Numeric. Optional override for on-target score cutoff.
#'        If NULL (default), auto-selects based on scoring method:
#'        - DeepCpf1/DeepSpCas9: 50 (percentage scale)
#'        - RuleSet1/Azimuth/enPAM+GB/DeepHF: 0.5 (probability-like scale)
#'        - RuleSet3: 0.1 (z-score scale)
#' @param offtarget Logical. Whether to run genome-wide off-target analysis
#'        (default TRUE). Set to FALSE for faster runs during prototyping.
#' @param ot_backend Character. Off-target search backend:
#'        \describe{
#'          \item{"bowtie"}{(Default) Fast mismatch-only search via crisprBowtie.
#'                Requires a one-time Bowtie index build. No bulge detection.}
#'          \item{"hybrid"}{Bowtie mismatch search on all gRNAs (<1 min),
#'                followed by CRISPRitz bulge refinement on recommended-pair gRNAs
#'                only. Best balance of speed and sensitivity. Falls back gracefully
#'                to Bowtie-only if CRISPRitz is unavailable. NOTE: This is an experimental
#'                feature - while the hybrid scoring approach works, it generally leads to
#'                no recommended pairs being returned due to stringency of filtering.}
#'          \item{"crispritz"}{Single-pass indexed CRISPRitz (bMax=1) on all
#'                gRNAs. No Bowtie dependency. Includes mismatch + bulge
#'                off-targets. Slower than hybrid (~100 min for TP53-scale)
#'                but fully self-contained. Use \code{full_bulge_scan = TRUE}
#'                for exhaustive bMax>=2 analysis (very slow). NOTE: *NOT RECOMMENDED* for
#'                local use. Still highly experimental, enormously memory-intensive, and not
#'                fully architected.}
#'        }
#'        Only used when \code{offtarget = TRUE}.
#' @param ot_threads Integer. Number of threads for CRISPRitz off-target searching
#'        (default 4). Ignored by Bowtie backend.
#' @param ot_timeout Integer. Maximum seconds for CRISPRitz execution (default 600).
#'        Ignored by Bowtie backend.
#' @param max_mismatches Integer. Maximum mismatches for off-target search (default 3).
#'        Clamped to 3 for Bowtie backend. CRISPRitz supports up to 6.
#' @param max_dna_bulges Integer. Maximum DNA bulge size (default 2).
#'        Used by CRISPRitz/hybrid backends only.
#' @param max_rna_bulges Integer. Maximum RNA bulge size (default 2).
#'        Used by CRISPRitz/hybrid backends only.
#' @param ot_canonical_only Logical (default TRUE). If TRUE, for Bowtie/hybrid backends
#'        search only canonical PAM off-targets. Matches CRISPRitz behaviour and MIT specificity
#'        calibration. Set to FALSE for analysis with non-canonical PAMs.
#' @param ot_detail_level Character. Controls retention of off-target hit details
#'        (default \code{"compact"}). Passed to \code{\link{score_offtargets}}.
#'        \describe{
#'          \item{"compact"}{(Default) Retains significant hits only (CFD > 0.01,
#'                top 200 per gRNA). Good balance of information and memory.}
#'          \item{"full"}{Retains all hits. WARNING: may exceed 1 GB for large genomes.}
#'          \item{"file"}{Writes full details to compressed RDS in cache directory;
#'                attaches compact version in memory. Path available via
#'                \code{result$offtarget_details_path}.}
#'          \item{"none"}{Discards hit details. Smallest footprint; per-gRNA
#'                specificity scores are still computed.}
#'        }
#' @param full_bulge_scan Logical (default FALSE). When TRUE and
#'        \code{ot_backend = "crispritz"}, passes the user's full bulge
#'        parameters (\code{max_dna_bulges}, \code{max_rna_bulges}) to
#'        CRISPRitz instead of capping at bMax=1. This is computationally
#'        expensive via WSL and unstable at bMax>=2 for multi-guide batches
#'        (sequential single-guide fallback will be triggered; ~18 min/guide).
#'        Ignored for \code{"bowtie"} and \code{"hybrid"} backends.
#' @param min_pair_specificity Numeric. Minimum pair-level specificity score
#'        (0-100, default 10) for recommended status, using the additive off-target
#'        burden model. See \code{\link{assemble_grna_pairs}} for details. Only
#'        applied when \code{offtarget = TRUE}.
#' @param design_primers Logical. Whether to design genotyping primers (default TRUE).
#' @param primer_max_wt Integer. Max WT amplicon size before switching to dual-pair strategy (default 3000).
#' @param primer_tm Numeric. Target melting temperature for primers (default 60.0).
#' @param top_n Integer. Number of top recommended pairs to plot (default 10; NULL = all).
#' @param quiet Logical. Suppress intermediate messages (default FALSE).
#' @param plot_mode Character. One of "heat" (default) or "arc". Passed to \code{plot_grna_design()}
#' @param interactive Logical; default FALSE. Use interactive plotly-based viewer.
#'
#' @return A named list with elements:
#'   - gene_id
#'   - gene_symbol
#'   - transcript_id
#'   - nuclease
#'   - exons (GRanges)
#'   - scored_grnas (GRanges)
#'   - pairs (data.frame with gRNAs and primers)
#'   - offtarget_details (data.frame, filtered per ot_detail_level; NULL if offtarget=FALSE)
#'   - offtarget_details_path (character, filepath when ot_detail_level="file"; NULL otherwise)
#'   - ot_detail_level (character, the detail level used)
#'   - plot (ggplot or plotly object)
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Cas9 workflow (default scoring: RuleSet1)
#' TP53_cas9 <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas9"
#' )
#'
#' # Cas9 workflow with DeepHF scoring (U6 promoter context)
#' TP53_deephf <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas9",
#'   score_method = "deephf",
#'   deephf_var = "wt_u6"
#' )
#'
#' # Cas9 workflow with DeepHF for RNP delivery (T7/synthetic gRNA)
#' TP53_deephf_rnp <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas9",
#'   score_method = "deephf",
#'   deephf_var = "wt_t7"
#' )
#'
#' # Cas9 workflow with DeepHF for eSpCas9(1.1)
#' TP53_espcas9 <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas9",
#'   score_method = "deephf",
#'   deephf_var = "esp"
#' )
#'
#' # Wild-type Cas12a workflow
#' TP53_cas12a <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas12a"
#' )
#'
#' # Engineered enCas12a workflow (expanded PAM)
#' TP53_encas12a <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "enCas12a"
#' )
#' # --- Off-target backend examples ---
#'
#' # Default: fast Bowtie mismatch-only (recommended for most users)
#' TP53_bowtie <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas9",
#'   ot_backend = "bowtie"       # default
#' )
#'
#' # Hybrid: Bowtie for all guides, CRISPRitz bulge refinement for recommended pairs
#' TP53_hybrid <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas9",
#'   ot_backend = "hybrid",
#'   max_dna_bulges = 2L,
#'   max_rna_bulges = 2L
#' )
#'
#' # CRISPRitz only: comprehensive search with bulges on all guides (slow)
#' TP53_crispritz <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas9",
#'   ot_backend = "crispritz",
#'   max_mismatches = 4L,
#'   max_dna_bulges = 2L,
#'   max_rna_bulges = 2L,
#'   ot_threads = 8L,
#'   ot_timeout = 1200L
#' )
#'
#' # Skip off-target scoring entirely (fast prototyping)
#' TP53_fast <- run_mutateR(
#'   gene_id = "TP53",
#'   species = "hsapiens",
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease = "Cas9",
#'   offtarget = FALSE
#' )
#' }
#'
#' @export
run_mutateR <- function(gene_id,
                        species,
                        genome,
                        nuclease = "Cas9",
                        transcript_id = NULL,
                        score_method = NULL,
                        tracr = "Chen2013",
                        deephf_var = c("wt_u6", "wt_t7", "wt", "esp", "hf"),
                        min_score = NULL,
                        offtarget = TRUE,
                        ot_backend = c( "bowtie", "hybrid", "crispritz"),
                        ot_threads = 4L,
                        ot_timeout = 600L,
                        max_mismatches = 3L,
                        max_dna_bulges = 2L,
                        max_rna_bulges = 2L,
                        ot_canonical_only = TRUE,
                        ot_detail_level = c("compact", "full", "file", "none"),
                        full_bulge_scan = FALSE,
                        min_pair_specificity = 10,
                        design_primers = TRUE,
                        primer_max_wt = 3000,
                        primer_tm = 60.0,
                        top_n = 10,
                        quiet = FALSE,
                        plot_mode = c("heat", "arc"),
                        interactive = FALSE) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
  })

  plot_mode <- match.arg(plot_mode)
  deephf_var <- match.arg(deephf_var)
  ot_backend <- match.arg(ot_backend)
  ot_detail_level <- match.arg(ot_detail_level)

  # ---- Normalise nuclease argument (string or NucleaseSpec) ----
  if (is.character(nuclease)) {
    nuclease <- match.arg(nuclease, c("Cas9", "Cas12a", "enCas12a"))
    nuc_spec  <- resolve_nuclease(nuclease)
    nuc_label <- nuclease
    is_custom_spec <- FALSE
  } else if (inherits(nuclease, "NucleaseSpec")) {
    nuc_spec  <- nuclease
    nuc_label <- if (isTRUE(nuc_spec$is_canonical)) nuc_spec$canonical_key else nuc_spec$name
    is_custom_spec <- !isTRUE(nuc_spec$is_canonical)
    # Re-bind `nuclease` to the canonical string when the spec is canonical,
    # so the existing switch() calls below still dispatch correctly.
    if (!is_custom_spec) nuclease <- nuc_spec$canonical_key
  } else {
    stop("`nuclease` must be a character string ('Cas9', 'Cas12a', 'enCas12a') ",
         "or a NucleaseSpec object built with nuclease_spec().")
  }

  if (!inherits(genome, "BSgenome"))
    stop("Please supply a valid BSgenome object.")

  # ---- Custom-spec gating: force-skip scoring stages ----
  if (is_custom_spec) {
    if (!quiet) {
      message("Custom NucleaseSpec '", nuc_spec$name,
              "' detected. On-target and off-target scoring will be skipped; ",
              "for lightweight scanning use run_mutateR_custom().")
    }
    score_method <- "none"
    offtarget    <- FALSE
  }

  # ---- Set default scoring method based on nuclease ----
  if (is.null(score_method)) {
    score_method <- switch(nuclease,
                           "Cas9" = "ruleset1",
                           "Cas12a" = "deepcpf1",
                           "enCas12a" = "enpamgb"
    )
    if (!quiet) message("Auto-selected scoring method '", score_method, "' for ", nuclease)
  }

  # ---- Validate DeepHF variant parameter usage ----
  if (score_method != "deephf" && deephf_var != "wt_u6") {
    # User explicitly set deephf variant but isn't using DeepHF
    if (!quiet) message("Note: 'deephf_var' parameter is only used with score_method='deephf'. Ignoring.")
  }

  ## ----- Step 1: Gene info -----
  if (!quiet) message("Retrieving gene/transcript information...")
  tx_info <- get_gene_info(gene_id, species)

  gene_symbol <- NA_character_
  if (!is.null(tx_info$canonical) && "external_gene_name" %in% names(tx_info$canonical)) {
    gene_symbol <- unique(tx_info$canonical$external_gene_name)[1]
  } else if ("external_gene_name" %in% names(tx_info$all)) {
    gene_symbol <- unique(tx_info$all$external_gene_name)[1]
  }

  canonical_tx <- if (is.null(transcript_id)) {
    if (!is.null(tx_info$canonical))
      tx_info$canonical$ensembl_transcript_id[1]
    else tx_info$all$ensembl_transcript_id[1]
  } else transcript_id

  if (!quiet) message("Using transcript: ", canonical_tx, " for gene: ", gene_id)

  ## ----- Step 2: Exon structures -----
  exons_gr <- get_exon_structures(canonical_tx, species, output = "GRanges")

  generate_early_plot <- function(pairs = NULL) {
    if (interactive) {
      plot_grna_interactive(exons_gr, pairs, transcript_id = canonical_tx,
                            gene_symbol = gene_symbol, species = species)
    } else {
      plot_grna_design(exons_gr, pairs, transcript_id = canonical_tx,
                       gene_symbol = gene_symbol, species = species,
                       mode = plot_mode)
    }
  }

  ## ----- Step 3A: Find gRNAs -----
  if (!quiet) {
    pam_info <- if (is_custom_spec) nuc_spec$pam else switch(nuclease,
                       "Cas9" = "NGG",
                       "Cas12a" = "TTTV",
                       "enCas12a" = "TTTN"
    )
    message("Locating ", nuc_spec$name, " target sites (PAM: ", pam_info, ")...")
  }

  hits <- if (is_custom_spec) {
    find_custom_grna_sites(exons_gr, genome, nuc_spec, require_full_context = FALSE)
  } else {
    switch(nuclease,
           "Cas9"     = find_cas9_sites(exons_gr, genome),
           "Cas12a"   = find_cas12a_sites(exons_gr, genome, pam = "TTTV"),
           "enCas12a" = find_cas12a_sites(exons_gr, genome, pam = "TTTN"))
  }

  if (is.null(hits) || length(hits) == 0) {
    warning("No gRNA sites identified for ", gene_id, " / ", nuc_label)
    return(list(gene_id = gene_id,
                gene_symbol = gene_symbol,
                transcript_id = canonical_tx,
                nuclease = nuc_label,
                nuclease_spec = nuc_spec,
                exons = exons_gr,
                scored_grnas = NULL,
                pairs = data.frame(),
                plot = generate_early_plot(NULL)))
  }

  # --- Step 3B: PAM distribution reporting ---
  if (!quiet && "pam_sequence" %in% names(mcols(hits))) {
    pam_counts <- table(hits$pam_sequence)
    message("PAM distribution:")
    print(pam_counts)
  }

  ## ----- Step 4A: On-target scoring -----
  if (is_custom_spec) {
    # Custom NucleaseSpec: skip scoring; annotate columns with NA so
    # downstream consumers see a consistent schema.
    scored_grnas <- hits
    GenomicRanges::mcols(scored_grnas)$ontarget_score <- NA_real_
    GenomicRanges::mcols(scored_grnas)$scoring_method <- "none"
    GenomicRanges::mcols(scored_grnas)$gc <- vapply(
      as.character(GenomicRanges::mcols(scored_grnas)$sequence_context),
      function(s) {
        if (is.na(s)) return(NA_real_)
        chars <- unlist(strsplit(toupper(s), ""))
        mean(chars %in% c("G", "C"))
      },
      numeric(1), USE.NAMES = FALSE
    )
  } else {
    if (!quiet) {
      msg <- paste0("Scoring gRNAs using model: ", score_method)
      if (score_method == "ruleset3") msg <- paste0(msg, " (tracrRNA: ", tracr, ")")
      if (score_method == "deephf")   msg <- paste0(msg, " (variant: ", deephf_var, ")")
      message(msg)
    }

    scored_grnas <- if (quiet) {
      suppressMessages(score_grnas(hits, method = score_method,
                                   tracr = tracr, deephf_var = deephf_var))
    } else {
      score_grnas(hits, method = score_method,
                  tracr = tracr, deephf_var = deephf_var)
    }
  }

  ## ----- Step 4B: Off-target scoring -----
  if (offtarget) {

    # CRISPRitz backend strategy (default: full_bulge_scan = FALSE):
    # Use indexed bMax=1 search on all guides — faster than brute-force FASTA
    # and captures bulge off-targets, eliminating the need for a Phase 2 pass.
    # Hybrid/bowtie: pass user params unchanged.
    # full_bulge_scan = TRUE: pass full user params (very slow at bMax>=2).
    cap_bulges        <- (ot_backend == "crispritz" && !full_bulge_scan)
    phase1_dna_bulges <- if (cap_bulges) min(max_dna_bulges, 1L) else max_dna_bulges
    phase1_rna_bulges <- if (cap_bulges) min(max_rna_bulges, 1L) else max_rna_bulges
    phase1_max_bulge  <- max(phase1_dna_bulges, phase1_rna_bulges)

    # Adaptive timeout: empirical model (bMax=1, 4t, hg38) = 34s + 15s/guide × 2.5 safety
    phase1_timeout <- ot_timeout
    if (cap_bulges && phase1_max_bulge > 0L) {
      n_guides_est   <- length(scored_grnas)
      phase1_timeout <- as.integer(max(ot_timeout,
                                       ceiling((34 + 15 * n_guides_est) * 2.5)))
    }

    if (!quiet) {
      if (ot_backend == "crispritz" && full_bulge_scan) {
        message("Full bulge scan on all gRNAs (backend: crispritz)...")
        message("  WARNING: This may take a long time for genes with many gRNAs.")
        if (max(max_dna_bulges, max_rna_bulges) >= 2L) {
          message("  WARNING: bMax>=2 is unstable for multi-guide batches.")
          message("  Sequential fallback will trigger on failure.")
        }
      } else if (ot_backend == "crispritz" && phase1_max_bulge > 0L) {
        est_min <- round((34 + 15 * length(scored_grnas)) / 60, 1)
        message("Running indexed off-target analysis (bMax=1) for all gRNAs ",
                "(backend: crispritz)...")
        message("  Includes mismatch + bulge off-targets in a single pass.")
        message("  Estimated time: ~", est_min, " min (timeout: ", phase1_timeout, "s)")
      } else if (ot_backend == "crispritz") {
        message("Running mismatch-only off-target analysis (backend: crispritz)...")
        message("  Note: Brute-force FASTA search; may be slow for many guides.")
      } else {
        message("Running genome-wide off-target analysis (backend: ", ot_backend, ")...")
      }
    }

    scored_grnas <- tryCatch(
      score_offtargets(scored_grnas,
                       genome         = genome,
                       nuclease       = nuclease,
                       exon_gr        = exons_gr,
                       backend        = ot_backend,
                       max_mismatches = max_mismatches,
                       max_dna_bulges = phase1_dna_bulges,
                       max_rna_bulges = phase1_rna_bulges,
                       canonical_only = ot_canonical_only,
                       detail_level   = ot_detail_level,
                       threads        = ot_threads,
                       timeout        = phase1_timeout,
                       quiet          = quiet),
      error = function(e) {
        warning("Off-target scoring failed: ", e$message,
                "\nContinuing without off-target scores.")
        scored_grnas
      }
    )

    # ---- Step 4C: Ensure OT columns exist (defensive) ----
    # If Phase 1 failed, scored_grnas lacks specificity columns.
    # Initialise with NA so downstream code doesn't crash.
    if (!"specificity_score" %in% names(mcols(scored_grnas))) {
      if (!quiet) message("Note: Off-target scores unavailable. ",
                          "Initialising specificity columns with NA.")
      mcols(scored_grnas)$specificity_score    <- NA_real_
      mcols(scored_grnas)$n_offtargets         <- NA_integer_
      mcols(scored_grnas)$has_exonic_offtarget <- FALSE
    }
  }

  ## ----- Step 5: Assembly -----
  if (!quiet) message("Assembling valid gRNA pairs for ", gene_id, " ...")

  valid_grnas <- suppressMessages({
    suppressWarnings({
      tmp <- capture.output(
        val <- filter_valid_grnas(exons_gr, genome, species,
                                  nuclease = nuc_spec,
                                  score_method = score_method,
                                  scored_grnas = scored_grnas,
                                  tracr = tracr,
                                  deephf_var = deephf_var)
      )
      val
    })
  })

  pairs_df <- tryCatch(
    assemble_grna_pairs(valid_grnas,
                        exon_gr = exons_gr,
                        transcript_id = canonical_tx,
                        species = species,
                        score_cutoff = min_score,
                        min_pair_specificity = min_pair_specificity),
    error = function(e) {
      warning("Pair assembly failed: ", e$message)
      NULL
    }
  )

  ## ----- Step 6A: Intragenic mode detection -----
  if (isTRUE(attr(pairs_df, "intragenic_mode")) && !quiet) {
    message("Detected intragenic assembly mode (≤2 exons).")
  }

  ## ----- Step 6B: Bulge refinement (Phase 2) -----
  # Phase 2 applies to hybrid backend only: Bowtie Phase 1 provides
  # mismatch-only scores; CRISPRitz Phase 2 adds bulge refinement for
  # recommended-pair gRNAs.
  # CRISPRitz backend: Phase 1 already runs indexed bMax=1 on ALL guides,
  # so Phase 2 would re-search the same subset at the same bMax — redundant.
  # Bowtie backend: mismatch-only by design, no Phase 2.
  needs_phase2 <- ot_backend == "hybrid"

  if (offtarget && needs_phase2 &&
      !is.null(pairs_df) && is.data.frame(pairs_df) &&
      nrow(pairs_df) > 0 && sum(pairs_df$recommended, na.rm = TRUE) > 0) {

    # ---- Phase 2 bulge depth selection ----
    # bMax=2 indexed search is unstable for multi-guide batches and
    # ~22x slower per guide than bMax=1. Cap Phase 2 at bMax=1 and
    # use the linear scaling model for adaptive timeout.
    phase2_dna_bulges <- min(max_dna_bulges, 1L)
    phase2_rna_bulges <- min(max_rna_bulges, 1L)

    if (max(max_dna_bulges, max_rna_bulges) > 1L && !quiet) {
      message("Note: Phase 2 bulge search capped at bMax=1 for stability and speed.")
      message("  bMax=2 indexed search is unstable for multi-guide batches (~22x slower,")
      message("  crashes with >1 guide). Use full_bulge_scan=TRUE for exhaustive bMax=2")
      message("  analysis (runs single-guide sequential; very slow).")
    }

    # ---- Phase 2 adaptive timeout ----
    # Empirical scaling model (bMax=1, 4 threads, hg38):
    #   time = 34s fixed + 15s/guide
    # Apply 2.5x safety margin for genome size variation and system load.
    rec_pairs <- pairs_df[pairs_df$recommended == TRUE, , drop = FALSE]
    n_phase2_spacers <- length(unique(c(
      as.character(rec_pairs$protospacer_sequence_5p),
      as.character(rec_pairs$protospacer_sequence_3p)
    )))

    phase2_timeout <- as.integer(max(
      ot_timeout,
      ceiling((34 + 15 * n_phase2_spacers) * 2.5)
    ))

    if (!quiet) {
      est_minutes <- round((34 + 15 * n_phase2_spacers) / 60, 1)
      message("Phase 2: refining ", sum(rec_pairs$recommended, na.rm = TRUE),
              " recommended pairs (", n_phase2_spacers, " unique spacers) ",
              "with CRISPRitz bulge search...")
      message("  Bulge depth: bDNA=", phase2_dna_bulges,
              ", bRNA=", phase2_rna_bulges)
      message("  Estimated time: ~", est_minutes, " min ",
              "(timeout: ", phase2_timeout, "s)")
    }

    refinement <- tryCatch(
      refine_with_bulges(
        pairs_df        = pairs_df,
        scored_grnas    = scored_grnas,
        genome          = genome,
        nuclease        = nuclease,
        exon_gr         = exons_gr,
        max_mismatches  = max_mismatches,
        max_dna_bulges  = phase2_dna_bulges,
        max_rna_bulges  = phase2_rna_bulges,
        bulge_penalty   = 0.2,
        threads         = ot_threads,
        timeout         = phase2_timeout,
        min_pair_specificity = min_pair_specificity,
        score_cutoff    = min_score,
        cache_dir       = NULL,
        detail_level    = ot_detail_level,
        quiet           = quiet
      ),
      error = function(e) {
        warning("Phase 2 bulge refinement failed: ", e$message,
                "\nRetaining Phase 1 specificity scores.")
        if (!quiet) message("Phase 2 fallback: retaining Bowtie Phase 1 specificity scores. ",
                            "Reason: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(refinement)) {
      old_details      <- attr(scored_grnas, "offtarget_details")
      old_details_path <- attr(scored_grnas, "offtarget_details_path")
      pairs_df         <- refinement$pairs
      scored_grnas     <- refinement$scored_grnas
      if (is.null(attr(scored_grnas, "offtarget_details")))
        attr(scored_grnas, "offtarget_details") <- old_details
      if (is.null(attr(scored_grnas, "offtarget_details_path")))
        attr(scored_grnas, "offtarget_details_path") <- old_details_path
    }
  }

  ## ----- Step 7: Primer design -----
  if (design_primers && !is.null(pairs_df) && nrow(pairs_df) > 0) {

    # Check if we have any recommended pairs
    n_rec <- sum(pairs_df$recommended, na.rm = TRUE)

    if (n_rec > 0) {
      if (!quiet) message("Designing genotyping primers for ", n_rec, " recommended pairs...")

      # Split dataframe
      df_rec  <- pairs_df[pairs_df$recommended == TRUE, ]
      df_rest <- pairs_df[pairs_df$recommended == FALSE, ]

      # Run design ONLY on recommended
      df_rec <- tryCatch({
        get_genotyping_primers(
          df_rec,
          exons_gr,
          genome,
          max_wt_amplicon = primer_max_wt,
          target_tm = primer_tm
        )
      }, error = function(e) {
        warning("Primer design failed: ", e$message)
        df_rec # return un-modified if fail
      })

      # Recombine safely using dplyr::bind_rows (handles empty df_rest automatically)
      pairs_df <- dplyr::bind_rows(df_rec, df_rest)

      # Sort recommended to top
      pairs_df <- pairs_df[order(pairs_df$recommended, decreasing = TRUE), ]

    } else {
      if (!quiet) message("No recommended pairs found; skipping primer design.")
    }
  }

  ## ----- Step 8: Plotting -----
  plot_obj <- NULL
  try({
    plot_pairs <- if (is.null(pairs_df)) data.frame() else pairs_df

    if (nrow(plot_pairs) == 0) {
      warning("No gRNA pairs met scoring/compatibility criteria.")
    }

    if (interactive) {
      plot_obj <- plot_grna_interactive(exons_gr,
                                        plot_pairs,
                                        transcript_id = canonical_tx,
                                        gene_symbol = gene_symbol,
                                        species = species)
    } else {
      plot_obj <- plot_grna_design(exons_gr,
                                   plot_pairs,
                                   gene_symbol = gene_symbol,
                                   transcript_id = canonical_tx,
                                   species = species,
                                   top_n = top_n,
                                   mode = plot_mode)
    }
  })

  ## ----- Step 9: Return -----
  if (!quiet) {
    ot_summary <- if (!offtarget) {
      "skipped"
    } else if (ot_backend == "crispritz" && full_bulge_scan) {
      paste0("crispritz/full-bulge (bMax=", max(max_dna_bulges, max_rna_bulges), ")")
    } else if (ot_backend == "crispritz") {
      paste0("crispritz/indexed (bMax=",
             min(max(max_dna_bulges, max_rna_bulges), 1L), ")")
    } else if (ot_backend == "hybrid") {
      phase2_bmax <- min(max(max_dna_bulges, max_rna_bulges), 1L)
      paste0("hybrid/two-phase (Phase2: bMax=", phase2_bmax, ")")
    } else {
      ot_backend
    }

    message("mutateR pipeline completed for ", gene_id,
            " (", nuc_label, "/", score_method,
            if (score_method == "deephf") paste0("/", deephf_var) else "",
            ", off-target: ", ot_summary,
            "), finding ",
            ifelse(is.null(pairs_df), 0, nrow(pairs_df)),
            " gRNA pairs.")
  }

  return(list(
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    transcript_id = canonical_tx,
    nuclease = nuc_label,
    nuclease_spec = nuc_spec,
    score_method = score_method,
    deephf_var = if (score_method == "deephf") deephf_var else NA_character_,
    ot_backend = if (offtarget) ot_backend else NA_character_,
    exons = exons_gr,
    scored_grnas = scored_grnas,
    offtarget_details      = if (offtarget) attr(scored_grnas, "offtarget_details") else NULL,
    offtarget_details_path = if (offtarget) attr(scored_grnas, "offtarget_details_path") else NULL,
    ot_detail_level        = if (offtarget) ot_detail_level else NA_character_,
    pairs = pairs_df,
    plot = plot_obj
  ))
}
