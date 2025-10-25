#' Map protein domain annotations to exons
#'
#' Retrieves annotated protein domains (currently InterPro only - N.B. this will be changed) for a given
#' Ensembl transcript, including amino acid coordinates, and assigns them to the
#' exons of the transcript. Each overlapping domain annotation is returned once,
#' with duplicates removed.
#'
#' @param transcript_id Character. Ensembl transcript ID.
#' @param species Character. Ensembl species short name (e.g., "hsapiens").
#'
#' @return A data.frame with columns:
#'   - ensembl_exon_id
#'   - exon_rank
#'   - domain_source ("InterPro")
#'   - domain_id
#'   - domain_desc (InterPro annotation)
#'   - domain_start / domain_end (amino acid coordinates)
#'
#' @export
map_protein_domains <- function(transcript_id, species) {
  if (missing(transcript_id)) {
    stop("You must provide an Ensembl transcript ID.")
  }
  if (missing(species)) {
    stop("You must provide a species in short name format (e.g., 'hsapiens').")
  }
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("The 'biomaRt' package is required.")
  }
  
  ensembl <- biomaRt::useEnsembl(biomart = "genes",
                                 dataset = paste0(species, "_gene_ensembl"))
  
  # Only InterPro attributes
  domain_data <- biomaRt::getBM(
    attributes = c("ensembl_exon_id",
                   "rank",
                   "interpro",
                   "interpro_short_description",
                   "interpro_start",
                   "interpro_end"),
    filters = "ensembl_transcript_id",
    values = transcript_id,
    mart = ensembl
  )
  
  if (nrow(domain_data) == 0) {
    warning("No InterPro domain information found for transcript: ", transcript_id, " in species: ", species)
    return(NULL)
  }
  
  # ----- Build InterPro dataframe -----
  interpro_df <- data.frame(
    ensembl_exon_id = domain_data$ensembl_exon_id,
    exon_rank       = domain_data$rank,
    domain_source   = "InterPro",
    domain_id       = domain_data$interpro,
    domain_desc     = domain_data$interpro_short_description,
    domain_start    = domain_data$interpro_start,
    domain_end      = domain_data$interpro_end,
    stringsAsFactors = FALSE
  )
  
  # ----- Remove NA / blank domain IDs -----
  interpro_df <- interpro_df[!is.na(interpro_df$domain_id) & interpro_df$domain_id != "", ]
  
  # ----- Deduplicate rows -----
  interpro_df <- unique(interpro_df)
  
  # ----- Sort consistently -----
  interpro_df <- interpro_df[order(interpro_df$exon_rank, interpro_df$domain_start), ]
  rownames(interpro_df) <- NULL
  
  return(interpro_df)
}