#' Run the complete mutateR workflow
#'
#' Executes all major steps of the mutateR pipeline end‑to‑end:
#' 1. Gene/Transcript retrieval
#' 2. Exon phase mapping
#' 3. gRNA finding & scoring (Cas9 or Cas12a)
#' 4. Exon-deletion pair assembly (Auto-detects scoring thresholds unless overridden)
#' 5. Genotyping primer design (Runs only on recommended pairs)
#' 6. Visualisation
#'
#' @param gene_id Character. Gene symbol or Ensembl Gene ID (ENSG...).
#' @param species Character. e.g. "hsapiens", "mmusculus", "drerio".
#' @param genome  BSgenome object (e.g. BSgenome.Hsapiens.UCSC.hg38).
#' @param nuclease Character. One of "Cas9" or "Cas12a" (default "Cas9").
#' @param transcript_id Optional Ensembl transcript ID to override canonical.
#' @param score_method Character. On‑target scoring model:
#'        Cas9: "ruleset1", "azimuth", "deephf";  Cas12a: "deepcpf1".
#' @param min_score Numeric. Optional override for on-target score cutoff.
#'        If NULL (default), auto-selects (0.5 for Cas9 models, 50 for DeepCpf1).
#' @param design_primers Logical. Whether to design genotyping primers (default TRUE).
#' @param primer_max_wt Integer. Max WT amplicon size before switching to dual-pair strategy (default 3000).
#' @param primer_tm Numeric. Target melting temperature for primers (default 60.0).
#' @param top_n Integer. Number of top recommended pairs to plot (default 10; NULL = all).
#' @param quiet Logical. Suppress intermediate messages (default FALSE).
#' @param plot_mode Character. One of "heat" (default) or "arc". Passed to \code{plot_grna_design()}
#' @param interactive Logical; default FALSE. Use interactive plotly-based viewer.
#'
#' @return A named list with elements:
#'   - transcript_id
#'   - exons (GRanges)
#'   - scored_grnas (GRanges)
#'   - pairs (data.frame with gRNAs and primers)
#'   - plot (ggplot or plotly object)
#'
#' @export
run_mutateR <- function(gene_id,
                        species,
                        genome,
                        nuclease = c("Cas9", "Cas12a"),
                        transcript_id = NULL,
                        score_method = NULL,
                        min_score = NULL,
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
  nuclease <- match.arg(nuclease)

  if (!inherits(genome, "BSgenome"))
    stop("Please supply a valid BSgenome object.")

  if (is.null(score_method)) {
    score_method <- if (nuclease == "Cas9") "ruleset1" else "deepcpf1" # Need to amend now that DeepSpCas9 is implemented.
  }

  ## ----- Step 1: Gene Info -----
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

  ## ----- Step 3: Find gRNAs -----
  if (!quiet) message("Locating ", nuclease, " target sites...")
  if (nuclease == "Cas9")
    hits <- find_cas9_sites(exons_gr, genome)
  else
    hits <- find_cas12a_sites(exons_gr, genome)

  if (is.null(hits) || length(hits) == 0) {
    warning("No gRNA sites identified for ", gene_id, " / ", nuclease)
    return(list(gene_id = gene_id,
                gene_symbol = gene_symbol,
                transcript_id = canonical_tx,
                exons = exons_gr,
                scored_grnas = NULL,
                pairs = data.frame(),
                plot = generate_early_plot(NULL)))
  }

  ## ----- Step 4: Scoring -----
  if (!quiet) message("Scoring gRNAs using model: ", score_method)

  if (quiet) {
    suppressMessages({scored_grnas <- score_grnas(hits, method = score_method)})
  } else {
    scored_grnas <- score_grnas(hits, method = score_method)
  }

  ## ----- Step 5: Assembly -----
  if (!quiet) message("Assembling valid gRNA pairs for ", gene_id, " ...")

  valid_grnas <- suppressMessages({
    suppressWarnings({
      tmp <- capture.output(
        val <- filter_valid_grnas(exons_gr, genome, species, score_method)
      )
      val
    })
  })

  pairs_df <- tryCatch(
    assemble_grna_pairs(valid_grnas,
                        exon_gr = exons_gr,
                        transcript_id = canonical_tx,
                        species = species,
                        score_cutoff = min_score), # Passes user override or NULL
    error = function(e) {
      warning("Pair assembly failed: ", e$message)
      NULL
    }
  )

  ## ----- Step 6: Intragenic Detection & Unpacking -----
  intragenic_mode <- FALSE
  if (is.list(pairs_df) && "pairs" %in% names(pairs_df)) {
    pairs_df <- pairs_df$pairs
    intragenic_mode <- TRUE
    if (!quiet) message("Detected intragenic assembly mode (≤2 exons).")
  }

  ## ----- Step 7: Primer Design -----
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
  if (!quiet) message("mutateR pipeline completed for ", gene_id,
                      ", finding ",
                      ifelse(is.null(pairs_df), 0, nrow(pairs_df)),
                      " gRNA pairs.")
  return(list(
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    transcript_id = canonical_tx,
    exons = exons_gr,
    scored_grnas = scored_grnas,
    pairs = pairs_df,
    plot = plot_obj
  ))
}
