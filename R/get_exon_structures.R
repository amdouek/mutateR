#' Get exon structures for a specified transcript
#'
#' Retrieves exon coordinates, reading frame phases, transcript CDS length,
#' and exon-specific coding lengths for each exon of a transcript.
#'
#' @param transcript_id Character. Ensembl transcript ID.
#' @param species Character. Species short name (e.g., "hsapiens", "drerio", "mmusculus").
#' @param output Character. "data.frame" (default when this function is run in isolation) or "GRanges". Use the GRanges output for Cas effector PAM scanning.
#'
#' @return A data.frame or GRanges with:
#'   - ensembl_exon_id
#'   - chromosome_name
#'   - exon_chrom_start / exon_chrom_end
#'   - strand ("+","-","*")
#'   - rank
#'   - start_phase / end_phase
#'   - cds_start / cds_end (genomic positions of CDS within exon)
#'   - transcript_cds_length (total CDS length, repeated per exon)
#'   - exon_cds_length (calculated per exon, NA if non-coding, e.g. UTR)
#'
#' @examples
#' \dontrun{
#' tx_info <- get_gene_info("TP53", "hsapiens")
#' canonical_tx <- tx_info$canonical$ensembl_transcript_id[1]
#' exons_df <- get_exon_structures(canonical_tx, "hsapiens", output="data.frame")
#' exons_gr <- get_exon_structures(canonical_tx, "hsapiens", output="GRanges")
#' }
#' @export
get_exon_structures <- function(transcript_id, species, 
                                output = c("data.frame", "GRanges")) {
  if (missing(transcript_id) || missing(species)) {
    stop("You must provide both a transcript_id (see the `get_gene_info()` output) and species.")
  }
  output <- match.arg(output)
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("The 'biomaRt' package is required. Please install it via Bioconductor.")
  }
  
  ensembl <- biomaRt::useEnsembl(biomart = "genes",
                                 dataset = paste0(species, "_gene_ensembl"))
  
  exons <- biomaRt::getBM(
    attributes = c("ensembl_exon_id",
                   "chromosome_name",
                   "exon_chrom_start",
                   "exon_chrom_end",
                   "strand",
                   "rank",
                   "phase",
                   "end_phase",
                   "cds_start",
                   "cds_end",
                   "cds_length"),
    filters = "ensembl_transcript_id",
    values = transcript_id,
    mart = ensembl
  )
  
  if (nrow(exons) == 0) {
    warning("No exon data found for transcript ", transcript_id)
    return(NULL)
  }
  
  # ----- Rename columns -----
  names(exons)[names(exons) == "phase"] <- "start_phase"
  names(exons)[names(exons) == "cds_length"] <- "transcript_cds_length"
  
  # ----- Convert default numeric strand identifiers to regular +/-/* symbols -----
  exons$strand <- ifelse(exons$strand == 1, "+",
                         ifelse(exons$strand == -1, "-", "*"))
  exons$strand[is.na(exons$strand)] <- "*"
  
  # ----- Ensure columns have proper types -----
  exons$chromosome_name    <- as.character(exons$chromosome_name)
  exons$rank               <- as.integer(exons$rank)
  exons$start_phase        <- as.integer(exons$start_phase)
  exons$end_phase          <- as.integer(exons$end_phase)
  exons$cds_start          <- as.integer(exons$cds_start)
  exons$cds_end            <- as.integer(exons$cds_end)
  exons$transcript_cds_length <- as.integer(exons$transcript_cds_length)
  
  # ----- Compute exon-specific CDS length -----
  exons$exon_cds_length <- ifelse(
    is.na(exons$cds_start) | is.na(exons$cds_end),
    NA_integer_,
    as.integer(exons$cds_end - exons$cds_start + 1)
  )
  
  # ----- Sort rows by exon rank -----
  exons <- exons[order(exons$rank), ]
  
  # ----- Output as data.frame -----
  if (output == "data.frame") {
    return(exons)
  }
  
  # ----- Output as a GRanges object -----
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The 'GenomicRanges' package is required for GRanges output.")
  }
  exons <- exons[!is.na(exons$exon_chrom_start) & !is.na(exons$exon_chrom_end), ]
  
  gr <- GenomicRanges::GRanges(
    seqnames = exons$chromosome_name,
    ranges   = IRanges::IRanges(start = exons$exon_chrom_start,
                                end   = exons$exon_chrom_end),
    strand   = exons$strand
  )
  # ----- Add all metadata (except the basic GRanges fields) -----
  mcols(gr) <- exons[, setdiff(names(exons),
                               c("chromosome_name",
                                 "exon_chrom_start",
                                 "exon_chrom_end",
                                 "strand"))]
  
  return(gr)
}