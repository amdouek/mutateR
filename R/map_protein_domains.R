#' Map protein domain annotations to exons (Pfam or InterPro, Ensembl ≥115‑compatible)
#'
#' Retrieves annotated Pfam or InterPro domains for a given Ensembl transcript.
#' Domain amino‑acid coordinates are mapped to the specific exons they overlap,
#' yielding biologically faithful exon–domain assignments.
#'
#' @param transcript_id Character. Ensembl transcript ID.
#' @param species Character. Ensembl species short name (e.g., "hsapiens").
#' @param source Character. Either "interpro" (default) or "pfam".
#'
#' @return data.frame with exon‑specific domain annotations.
#' @export
#'
map_protein_domains <- function(transcript_id,
                                species,
                                source = c("interpro", "pfam")) {

  source <- tolower(source)
  source <- match.arg(source, c("interpro","pfam"))

  if (!requireNamespace("biomaRt", quietly = TRUE))
    stop("The 'biomaRt' package is required.")
  if (!requireNamespace("httr", quietly = TRUE))
    stop("The 'httr' package is required for Pfam name lookup via InterPro API.")

  ## ----- Helper to retrieve exon structures (for mapping) -----
  exons_df <- tryCatch(
    get_exon_structures(transcript_id, species, output = "data.frame"),
    error = function(e) NULL
  )
  if (is.null(exons_df) || nrow(exons_df) == 0)
    stop("Could not retrieve exon structures for ", transcript_id)

  # Build exon to amino acid position ranges
  exons_df <- exons_df[!is.na(exons_df$exon_cds_length), ]
  exons_df$cds_aa_len <- floor(exons_df$exon_cds_length / 3)
  exons_df$cum_start <- c(1, cumsum(head(exons_df$cds_aa_len, -1)) + 1)
  exons_df$cum_end   <- cumsum(exons_df$cds_aa_len)

  ## ----- Connect to Ensembl Genes mart -----
  ensembl <- biomaRt::useEnsembl(biomart = "genes",
                                 dataset = paste0(species, "_gene_ensembl"))
  safe_getBM <- function(attrs) {
    biomaRt::getBM(attributes = attrs,
                   filters = "ensembl_transcript_id",
                   values  = transcript_id,
                   mart    = ensembl)
  }
  attr_names <- suppressMessages(biomaRt::listAttributes(ensembl)$name)

  ## ---- Internal helper: InterPro REST fetch for Pfam names ----
  fetch_pfam_desc <- function(pfam_ids) {
    pfam_ids <- unique(pfam_ids[!is.na(pfam_ids) & pfam_ids != ""])
    if (length(pfam_ids) == 0) return(character())
    descs <- vapply(pfam_ids, function(pid) {
      url <- paste0("https://www.ebi.ac.uk/interpro/api/entry/pfam/", pid)
      resp <- tryCatch(httr::GET(url, httr::timeout(3)), error = function(e) NULL)
      if (is.null(resp) || httr::status_code(resp) != 200) return(NA_character_)
      parsed <- tryCatch(httr::content(resp, as = "parsed", type = "application/json"),
                         error = function(e) NULL)
      val <- if (!is.null(parsed$metadata$name)) parsed$metadata$name else NA_character_
      if (length(val) > 1) val <- val[1]
      as.character(val)
    }, character(1))
    stats::setNames(descs, pfam_ids)
  }

  ## ----- Domain retrieval (Pfam or InterPro) -----
  if (source == "pfam") {
    pfam_attrs <- c("ensembl_transcript_id","ensembl_peptide_id",
                    "pfam","pfam_start","pfam_end")
    valid_attrs <- pfam_attrs[pfam_attrs %in% attr_names]
    if (length(valid_attrs) < 3) {
      message("Pfam attributes unavailable; falling back to InterPro.")
      return(map_protein_domains(transcript_id, species, "interpro"))
    }
    message("Retrieving Pfam domain annotations from Ensembl Genes mart...")
    domain_data <- safe_getBM(valid_attrs)
    if (nrow(domain_data) == 0)
      return(map_protein_domains(transcript_id, species, "interpro"))

    # Friendly names
    pfam_map <- fetch_pfam_desc(unique(domain_data$pfam))
    domain_data$pfam_description <- ifelse(domain_data$pfam %in% names(pfam_map),
                                           pfam_map[domain_data$pfam],
                                           domain_data$pfam)
    out <- data.frame(
      domain_source = "pfam",
      domain_id   = domain_data$pfam,
      domain_desc = domain_data$pfam_description,
      domain_start = domain_data$pfam_start,
      domain_end   = domain_data$pfam_end,
      stringsAsFactors = FALSE
    )
  } else {
    interpro_attrs <- c("ensembl_transcript_id","ensembl_peptide_id",
                        "interpro","interpro_start","interpro_end",
                        "interpro_short_description")
    valid_attrs <- interpro_attrs[interpro_attrs %in% attr_names]
    message("Retrieving InterPro domain annotations from Ensembl Genes mart...")
    domain_data <- safe_getBM(valid_attrs)
    if (nrow(domain_data) == 0) return(NULL)
    out <- data.frame(
      domain_source = "interpro",
      domain_id   = domain_data$interpro,
      domain_desc = domain_data$interpro_short_description,
      domain_start = domain_data$interpro_start,
      domain_end   = domain_data$interpro_end,
      stringsAsFactors = FALSE
    )
  }

  ## ----- Map domain amino‑acid coordinates to exon ranks -----
  map_domain_to_exons <- function(dstart, dend) {
    hits <- exons_df$rank[
      (exons_df$cum_start <= dend) & (exons_df$cum_end >= dstart)
    ]
    if (length(hits) == 0) hits <- NA_integer_
    hits
  }

  expanded_rows <- lapply(seq_len(nrow(out)), function(i) {
    exons_hit <- map_domain_to_exons(out$domain_start[i], out$domain_end[i])
    if (is.na(exons_hit[1])) return(NULL)
    data.frame(
      exon_rank       = exons_hit,
      domain_source   = out$domain_source[i],
      domain_id       = out$domain_id[i],
      domain_desc     = out$domain_desc[i],
      domain_start    = out$domain_start[i],
      domain_end      = out$domain_end[i],
      stringsAsFactors = FALSE
    )
  })
  out_exp <- do.call(rbind, expanded_rows)

  # Add exon IDs from exon_df
  out_exp <- merge(out_exp,
                   exons_df[, c("ensembl_exon_id","rank")],
                   by.x = "exon_rank", by.y = "rank", all.x = TRUE)

  out_exp <- out_exp[, c("ensembl_exon_id","exon_rank","domain_source",
                         "domain_id","domain_desc","domain_start","domain_end")]
  out_exp <- unique(out_exp)
  rownames(out_exp) <- NULL
  return(out_exp)
}
