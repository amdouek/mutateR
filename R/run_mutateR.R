#' Run the complete mutateR workflow
#'
#' Executes all major steps of the mutateR pipeline end‑to‑end
#'
#' @param gene_id Character. Gene symbol or Ensembl Gene ID (ENSG...).
#' @param species Character. e.g. "hsapiens", "mmusculus", "drerio".
#' @param genome  BSgenome object (e.g. BSgenome.Hsapiens.UCSC.hg38).
#' @param nuclease Character. One of "Cas9" or "Cas12a" (default "Cas9").
#' @param transcript_id Optional Ensembl transcript ID to override canonical.
#' @param score_method Character. On‑target scoring model:
#'        Cas9: "ruleset1", "azimuth", "deephf";  Cas12a: "deepcpf1".
#' @param top_n Integer. Number of top recommended pairs to plot (default 10; NULL = all).
#' @param quiet Logical. Suppress intermediate messages (default FALSE).
#'
#' @return A named list with elements:
#'   - transcript_id
#'   - exons
#'   - scored_grnas
#'   - pairs
#'   - plot
#'
#' @export
run_mutateR <- function(gene_id,
                        species,
                        genome,
                        nuclease = c("Cas9", "Cas12a"),
                        transcript_id = NULL,
                        score_method = NULL,
                        top_n = 10,
                        quiet = FALSE) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
  })
  
  nuclease <- match.arg(nuclease)
  if (!inherits(genome, "BSgenome"))
    stop("Please supply a valid BSgenome object for the target species. \n Available genomes can be seen using BSgenome::available.genomes()")
  
  ## ----- Determine scoring model defaults -----
  if (is.null(score_method)) {
    score_method <- if (nuclease == "Cas9") "ruleset1" else "deepcpf1"
  }
  
  ## ----- Step 1: Retrieve canonical transcript -----
  if (!quiet) message("Retrieving gene/transcript information...")
  tx_info <- get_gene_info(gene_id, species)
  gene_symbol <- NA_character_
  if (!is.null(tx_info$canonical) &&
      "external_gene_name" %in% names(tx_info$canonical)) {
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
  
  ## ----- Step 3: Find gRNAs -----
  if (!quiet) message("Locating ", nuclease, " target sites...")
  if (nuclease == "Cas9")
    hits <- find_cas9_sites(exons_gr, genome)
  else
    hits <- find_cas12a_sites(exons_gr, genome)
  
  if (is.null(hits) || length(hits) == 0) {
    warning("No gRNA sites identified for ", gene_id, " / ", nuclease)
    return(list(transcript_id = canonical_tx,
                exons = exons_gr,
                scored_grnas = NULL,
                pairs = data.frame(),
                plot = plot_grna_design(exons_gr, NULL,
                                        transcript_id = canonical_tx)))
  }
  
  ## ----- Step 4: Calculate on-target gRNA scores -----
  if (!quiet) message("Scoring gRNAs using model: ", score_method)
  os_is_windows <- identical(.Platform$OS.type, "windows")
  if (os_is_windows && tolower(score_method) == "deepcpf1") {
    warning("DeepCpf1 model unsupported on Windows; skipping scoring.")
    scored_grnas <- hits
    mcols(scored_grnas)$gc <- NA_real_
    mcols(scored_grnas)$ontarget_score <- NA_real_
  } else {
    if (quiet) {
      suppressMessages({
        scored_grnas <- score_grnas(hits, method = score_method)
      })
    } else {
      scored_grnas <- score_grnas(hits, method = score_method)
    }
  }
  
  ## ----- Step 5: Filter and assemble gRNA pairs -----
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
                        species = species),
    error = function(e) {
      warning("Pair assembly failed: ", e$message)
      NULL
    }
  )
  
  ## ----- Step 6: Detect intragenic mode & normalise output (for genes with >=2 exons) -----
  intragenic_mode <- FALSE
  if (is.list(pairs_df) && "pairs" %in% names(pairs_df)) {
    pairs_df <- pairs_df$pairs
    intragenic_mode <- TRUE
    message("Detected intragenic assembly mode (≤2 exons).")
  }
  
  ## ----- Step 7: Plot generation -----
  if (is.null(pairs_df) ||
      (is.data.frame(pairs_df) && nrow(pairs_df) == 0)) {
    warning("No gRNA pairs met scoring/compatibility criteria.")
    plot_obj <- tryCatch({
      plot_grna_design(exons_gr,
                       if (intragenic_mode) data.frame() else NULL,
                       gene_symbol = gene_symbol,
                       transcript_id = canonical_tx,
                       top_n = NULL)
    }, error = function(e) {
      warning("Plot generation failed: ", e$message)
      NULL
    })
  } else {
    plot_obj <- tryCatch({
      plot_grna_design(exons_gr, pairs_df,
                       gene_symbol = gene_symbol,
                       transcript_id = canonical_tx,
                       top_n = top_n)
    }, error = function(e) {
      warning("Plot generation failed: ", e$message)
      NULL
    })
  }
  
  ## ----- Step 8: Wrap and return -----
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