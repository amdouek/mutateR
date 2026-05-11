#' @title Lightweight pipeline for custom nucleases (no scoring)
#'
#' @description Streamlined alternative to \code{\link{run_mutateR}} for users
#' who only need site-scanning and (optionally) deletion-pair assembly and
#' primer design with a user-defined nuclease. The pipeline runs:
#' \enumerate{
#'   \item Gene/transcript retrieval
#'   \item Exon-phase mapping
#'   \item Custom site finding via \code{\link{find_custom_grna_sites}}
#'   \item Optional pair assembly (phase-compatibility + deletion-size filters)
#'   \item Optional genotyping primer design
#'   \item Optional visualisation
#' }
#'
#' On-target and off-target scoring are not performed — those steps depend on
#' nuclease-specific trained models that are not available for arbitrary
#' effectors. If you want scoring on the canonical Cas9 / Cas12a / enCas12a
#' workflows, use \code{\link{run_mutateR}}.
#'
#' \code{run_mutateR()} also accepts a custom \code{NucleaseSpec} as its
#' \code{nuclease} argument — that variant follows the full assembly pipeline
#' (still skipping scoring). \code{run_mutateR_custom()} is the explicit
#' "scanning + assembly only" entry point with a smaller surface area and no
#' scoring/off-target parameters.
#'
#' @param gene_id Character. Gene symbol or Ensembl Gene ID.
#' @param species Character. e.g. "hsapiens".
#' @param genome  BSgenome object.
#' @param nuclease_spec A \code{NucleaseSpec} object built with
#'        \code{\link{nuclease_spec}}. Canonical specs are accepted but will
#'        \emph{not} trigger scoring on this pipeline — use
#'        \code{\link{run_mutateR}} if you want scoring.
#' @param transcript_id Optional Ensembl transcript ID to override canonical.
#' @param assemble_pairs Logical. If TRUE (default), run phase filtering and
#'        deletion-pair assembly. If FALSE, return only the per-site GRanges
#'        and skip pair assembly entirely.
#' @param design_primers Logical. Whether to design genotyping primers for
#'        assembled pairs (default FALSE; turn on if you want amplicon design).
#' @param primer_max_wt Integer. Max WT amplicon size before switching to
#'        dual-pair primer strategy (default 3000). Ignored unless
#'        \code{design_primers = TRUE}.
#' @param primer_tm Numeric. Target melting temperature for primers
#'        (default 60.0). Ignored unless \code{design_primers = TRUE}.
#' @param plot Logical. If TRUE (default), produce a ggplot visualisation of
#'        the assembled design.
#' @param plot_mode Character. One of "heat" (default) or "arc". Passed to
#'        \code{\link{plot_grna_design}}.
#' @param top_n Integer. Number of pairs to plot (default 10; NULL = all).
#' @param interactive Logical (default FALSE). Use plotly-based viewer.
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{gene_id}, \code{gene_symbol}, \code{transcript_id}}{Identifiers used.}
#'     \item{\code{nuclease}}{Character name of the supplied spec.}
#'     \item{\code{nuclease_spec}}{The resolved \code{NucleaseSpec} object.}
#'     \item{\code{exons}}{Exon \code{GRanges}.}
#'     \item{\code{sites}}{\code{GRanges} of candidate sites (output of
#'           \code{find_custom_grna_sites()}).}
#'     \item{\code{pairs}}{Data.frame of assembled deletion pairs
#'           (\code{NULL} when \code{assemble_pairs = FALSE}).}
#'     \item{\code{plot}}{ggplot/plotly object or \code{NULL}.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' aac <- nuclease_spec(
#'   name               = "AacCas12b",
#'   pam                = "TTN",
#'   pam_side           = "5prime",
#'   protospacer_length = 20L,
#'   cut_offset_top     = -5L,
#'   cut_offset_bottom  = -9L,
#'   activity           = "cut",
#'   grna_architecture  = "crRNA + tracr"
#' )
#'
#' result <- run_mutateR_custom(
#'   gene_id       = "TP53",
#'   species       = "hsapiens",
#'   genome        = BSgenome.Hsapiens.UCSC.hg38,
#'   nuclease_spec = aac
#' )
#' result$sites
#' result$pairs
#' }
#'
#' @export
run_mutateR_custom <- function(gene_id,
                               species,
                               genome,
                               nuclease_spec,
                               transcript_id  = NULL,
                               assemble_pairs = TRUE,
                               design_primers = FALSE,
                               primer_max_wt  = 3000,
                               primer_tm      = 60.0,
                               plot           = TRUE,
                               plot_mode      = c("heat", "arc"),
                               top_n          = 10,
                               interactive    = FALSE,
                               quiet          = FALSE) {

  plot_mode <- match.arg(plot_mode)

  spec <- resolve_nuclease(nuclease_spec)
  nuc_label <- if (isTRUE(spec$is_canonical)) spec$canonical_key else spec$name

  if (!inherits(genome, "BSgenome"))
    stop("Please supply a valid BSgenome object.")

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
    else
      tx_info$all$ensembl_transcript_id[1]
  } else {
    transcript_id
  }
  if (!quiet) message("Using transcript: ", canonical_tx, " for gene: ", gene_id)

  ## ----- Step 2: Exon structures -----
  exons_gr <- get_exon_structures(canonical_tx, species, output = "GRanges")

  ## ----- Step 3: Site finding (custom spec) -----
  if (!quiet) {
    message("Locating ", spec$name, " target sites (PAM: ", spec$pam,
            ", ", spec$pam_side, ", spacer ", spec$protospacer_length, " nt)...")
  }
  sites <- find_custom_grna_sites(exons_gr, genome, spec,
                                  require_full_context = FALSE)

  if (is.null(sites) || length(sites) == 0) {
    warning("No ", spec$name, " sites identified in ", gene_id, ".")
    return(list(
      gene_id        = gene_id,
      gene_symbol    = gene_symbol,
      transcript_id  = canonical_tx,
      nuclease       = nuc_label,
      nuclease_spec  = spec,
      exons          = exons_gr,
      sites          = NULL,
      pairs          = NULL,
      plot           = NULL
    ))
  }

  # Annotate scoring-skip placeholders so downstream functions see a
  # consistent schema. These are intentionally NA — no scoring was performed.
  GenomicRanges::mcols(sites)$ontarget_score <- NA_real_
  GenomicRanges::mcols(sites)$scoring_method <- "none"

  ## ----- Step 4: Optional pair assembly -----
  pairs_df <- NULL
  if (assemble_pairs) {
    if (!quiet) message("Assembling phase-compatible deletion pairs...")

    valid_grnas <- suppressMessages(
      filter_valid_grnas(exons_gr, genome, species,
                         nuclease     = spec,
                         score_method = "none",
                         scored_grnas = sites)
    )

    if (!is.null(valid_grnas) && length(valid_grnas) > 0) {
      pairs_df <- tryCatch(
        assemble_grna_pairs(valid_grnas,
                            exon_gr       = exons_gr,
                            transcript_id = canonical_tx,
                            species       = species),
        error = function(e) {
          warning("Pair assembly failed: ", e$message)
          NULL
        }
      )
    }

    ## ----- Step 5: Optional primer design -----
    if (design_primers &&
        !is.null(pairs_df) && nrow(pairs_df) > 0 &&
        sum(pairs_df$recommended, na.rm = TRUE) > 0) {

      if (!quiet) message("Designing genotyping primers for recommended pairs...")
      df_rec  <- pairs_df[pairs_df$recommended == TRUE, , drop = FALSE]
      df_rest <- pairs_df[pairs_df$recommended == FALSE, , drop = FALSE]

      df_rec <- tryCatch({
        get_genotyping_primers(df_rec, exons_gr, genome,
                               max_wt_amplicon = primer_max_wt,
                               target_tm       = primer_tm)
      }, error = function(e) {
        warning("Primer design failed: ", e$message)
        df_rec
      })

      pairs_df <- dplyr::bind_rows(df_rec, df_rest)
      pairs_df <- pairs_df[order(pairs_df$recommended, decreasing = TRUE), ]
    }
  }

  ## ----- Step 6: Plot -----
  plot_obj <- NULL
  if (plot) {
    try({
      plot_pairs <- if (is.null(pairs_df)) data.frame() else pairs_df
      if (interactive) {
        plot_obj <- plot_grna_interactive(exons_gr, plot_pairs,
                                          transcript_id = canonical_tx,
                                          gene_symbol   = gene_symbol,
                                          species       = species)
      } else {
        plot_obj <- plot_grna_design(exons_gr, plot_pairs,
                                     gene_symbol   = gene_symbol,
                                     transcript_id = canonical_tx,
                                     species       = species,
                                     top_n         = top_n,
                                     mode          = plot_mode)
      }
    })
  }

  if (!quiet) {
    n_sites <- length(sites)
    n_pairs <- if (is.null(pairs_df)) 0L else nrow(pairs_df)
    message("run_mutateR_custom completed: ", n_sites, " site(s), ",
            n_pairs, " candidate pair(s).")
  }

  list(
    gene_id        = gene_id,
    gene_symbol    = gene_symbol,
    transcript_id  = canonical_tx,
    nuclease       = nuc_label,
    nuclease_spec  = spec,
    exons          = exons_gr,
    sites          = sites,
    pairs          = pairs_df,
    plot           = plot_obj
  )
}
