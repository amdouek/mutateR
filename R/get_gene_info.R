#' @title Retrieve transcript information for a selected gene from a given species
#'
#' @description Queries Ensembl via biomaRt to retrieve transcript annotation for a given
#' gene symbol or Ensembl gene ID in a specified species. Canonical
#' transcript(s), if any, are flagged and the output separated into a list of
#' two data.frames (all transcripts, and just the canonical transcript
#' (if applicable)).
#'
#' @param gene_id Character. Identifier for selected gene of interest: either gene symbol
#'        (via external_gene_name) or an Ensembl gene ID (e.g., ENSG...).
#' @param species Character. Species short name (e.g., "hsapiens", "drerio", "mmusculus", etc.).
#' @param id_type Character. One of "symbol" (default) or "ensembl_gene_id".
#'
#' @return A list with:
#'   - canonical: data.frame of canonical transcript(s), or NULL if no transcript annotated as canonical for the selected gene
#'   - all: full data.frame of transcripts, canonical (if present) raised to top
#'
#' @examples
#' \dontrun{
#' tx_info <- get_gene_info("TP53", "hsapiens")
#' tx_info$canonical    # canonical transcript(s)
#' tx_info$all          # all transcripts
#' }
#' @export
get_gene_info <- function(gene_id, species, id_type = c("symbol", "ensembl_gene_id")) {
  id_type <- match.arg(id_type)

  if (missing(species)) {
    stop("You must specify a species using the short name format (e.g., 'hsapiens', 'drerio', 'mmusculus', etc.).")
  }
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("The 'biomaRt' package is required but not installed. Please install it via Bioconductor.")
  }

  ensembl <- biomaRt::useEnsembl(biomart = "genes",
                                 dataset = paste0(species, "_gene_ensembl")) # Connects to Ensembl

  filter <- if (id_type == "symbol") "external_gene_name" else "ensembl_gene_id"

  transcripts <- biomaRt::getBM(
    attributes = c("ensembl_gene_id",
                   "ensembl_transcript_id",
                   "chromosome_name",
                   "start_position",
                   "end_position",
                   "strand",
                   "external_gene_name",
                   "transcript_is_canonical"),
    filters = filter,
    values  = gene_id,
    mart    = ensembl
  )  # Retrieves specified information from Ensembl via biomaRt

  # ----- Convert canonical flag to logical -----
  transcripts$transcript_is_canonical <- transcripts$transcript_is_canonical == 1

  # ----- Convert numeric strand to symbols -----
  transcripts$strand <- ifelse(transcripts$strand == 1,
                               "+",
                               ifelse(transcripts$strand == -1, "-", "*"))

  # ----- Separate canonical transcript(s) -----
  canonical_df <- subset(transcripts, transcript_is_canonical == TRUE)

  # ----- Sort full table so canonical transcript, if present, is at top -----
  transcripts_sorted <- transcripts[order(!transcripts$transcript_is_canonical,
                                          transcripts$ensembl_transcript_id), ]

  # ----- Return both 'canonical' and 'all' data.frames as a list -----
  return(list(canonical = if (nrow(canonical_df) > 0) canonical_df else NULL,
              all = transcripts_sorted))
}
