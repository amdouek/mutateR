#' Integrate gRNA scoring data with exon‑pair compatibility
#'
#' Returns scored Cas9 gRNAs that lie within exons
#' participating in phase‑compatible exon pairs or predicted
#' terminal‑exon frameshift cases (permitted PTC deletions).
#'
#' @param exon_gr  GRanges from get_exon_structures(output="GRanges").
#' @param genome   BSgenome object for the relevant species.
#' @param species  Character, e.g. "hsapiens".
#' @param score_method Character. One of "ruleset1", "azimuth", or "deephf".
#'
#' @return GRanges of scored, biologically admissible gRNAs.
#' @export
filter_valid_grnas <- function(exon_gr,
                               genome,
                               species,
                               score_method = c("ruleset1","azimuth","deephf","deepcpf1","deepspcas9")) {

  score_method <- match.arg(score_method)

  ## ---- 0. Early‑exit for single‑exon and two‑exon transcripts ----------
  n_exons <- length(exon_gr)
  if (n_exons <= 2) {
    warning("Transcript has ", n_exons, " exon(s). Skipping phase‑compatibility checks.")

    # Branch based on method to determine nuclease
    if (score_method == "deepcpf1") {
      grna_sites <- find_cas12a_sites(exon_gr, genome)
    } else {
      grna_sites <- find_cas9_sites(exon_gr, genome)
    }

    if (is.null(grna_sites) || length(grna_sites) == 0) {
      warning("No target sites found.")
      return(NULL)
    }
    grna_sites <- score_grnas(grna_sites, method = score_method)
    return(grna_sites)
  }

  ## ---- 1. Exon info and pair classification --------------------------
  exon_meta <- as.data.frame(mcols(exon_gr))
  exon_meta$rank <- seq_len(nrow(exon_meta))

  comp_pairs <- check_exon_phase(exon_meta, include_contiguous = FALSE)

  # Assess frameshift/PTC status for each exon pair
  fs_list <- lapply(seq_len(nrow(comp_pairs)), function(i) {
    with(comp_pairs[i, ],
         check_frameshift_ptc(exon_meta,
                              exon_5p = exon_5p,
                              exon_3p = exon_3p))
  })
  fs_df <- do.call(rbind, lapply(fs_list, as.data.frame))
  pair_info <- cbind(comp_pairs, fs_df)

  # Accepted exons = appear in compatible pairs
  # OR are part of terminal‑exon PTC‑tolerated cases
  accepted_exons <- unique(c(
    pair_info$exon_5p[pair_info$compatible | pair_info$terminal_exon_case],
    pair_info$exon_3p[pair_info$compatible | pair_info$terminal_exon_case]
  ))

  ## ---- 2. Collect and score Cas9 sites -------------------------------
  if (score_method == "deepcpf1") {
    grna_sites <- find_cas12a_sites(exon_gr, genome)
  } else {
    grna_sites <- find_cas9_sites(exon_gr, genome)
  }

  if (is.null(grna_sites) || length(grna_sites) == 0) {
    warning("No sites found within provided exons.")
    return(NULL)
  }

  grna_sites <- score_grnas(grna_sites, method = score_method)

  ## ---- 3. Annotate exon rank if missing -------------------------------
  if (!"exon_rank" %in% names(mcols(grna_sites))) {
    hits <- GenomicRanges::findOverlaps(grna_sites, exon_gr)
    grna_sites$exon_rank <- NA_integer_
    grna_sites$exon_rank[queryHits(hits)] <- exon_meta$rank[subjectHits(hits)]
  }

  ## ---- 4. Filter to allowed exon ranks -------------------------------
  valid <- grna_sites[grna_sites$exon_rank %in% accepted_exons]

  message("Retained ", length(valid),
          " scored gRNAs within phase‑compatible or terminal‑exon‑PTC cases.")

  valid
}
