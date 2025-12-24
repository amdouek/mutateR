#' @title Integrate gRNA scoring data with exon‑pair compatibility
#'
#' @description Returns scored Cas9/Cas12a gRNAs that lie within exons
#' participating in phase‑compatible exon pairs or predicted
#' terminal‑exon PTC cases.
#'
#' @param exon_gr  GRanges. Returned by get_exon_structures(output="GRanges").
#' @param genome   BSgenome object for the relevant species.
#' @param species  Character, using short-form species name e.g. "hsapiens".
#' @param nuclease Character. One of "Cas9", "Cas12a" or "enCas12a".
#' @param score_method Character. One of "ruleset1", "azimuth", "deepspcas9",
#'        "deephf", "ruleset3" (Cas9), "deepcpf1" (Cas12a), or "enpamgb" (enCas12a).
#' @param scored_grnas Precomputed on-target scores passed from elsewhere (default NULL, evokes scoring).
#' @param tracr Character. Either "Chen2013" (default) or "Hsu2013" - only for ruleset3 scoring.
#' @param deephf_var Character. For DeepHF scoring; one of "wt", "wt_u6" (default), "wt_t7", "esp", or "hf".
#'        See \code{\link{recommend_deephf_model}} for guidance on model selection.
#'
#' @return GRanges of scored, mutateR-allowed gRNAs.
#' @export
filter_valid_grnas <- function(exon_gr,
                               genome,
                               species,
                               nuclease = c("Cas9", "Cas12a", "enCas12a"),
                               score_method = c("ruleset1", "azimuth", "deepspcas9",
                                                "deephf", "ruleset3",
                                                "deepcpf1", "enpamgb"),
                               scored_grnas = NULL,
                               tracr = "Chen2013",
                               deephf_var = c("wt_u6", "wt_t7", "wt", "esp", "hf")) {

  nuclease <- match.arg(nuclease)
  score_method <- match.arg(score_method)
  deephf_var <- match.arg(deephf_var)

  # ---- Validate nuclease/score_method compatibility ----
  cas9_methods <- c("ruleset1", "azimuth", "deepspcas9", "deephf", "ruleset3")
  cas12a_methods <- c("deepcpf1")
  encas12a_methods <- c("enpamgb")

  if (nuclease == "Cas9" && !score_method %in% cas9_methods) {
    warning("Score method '", score_method, "' is not designed for Cas9. ",
            "Consider using one of: ", paste(cas9_methods, collapse = ", "))
  } else if (nuclease == "Cas12a" && !score_method %in% c(cas12a_methods, encas12a_methods)) {
    warning("Score method '", score_method, "' is not designed for Cas12a. ",
            "Consider using one of: ", paste(cas12a_methods, collapse = ", "))
  } else if (nuclease == "enCas12a" && !score_method %in% encas12a_methods) {
    warning("Score method '", score_method, "' is not optimised for enCas12a. ",
            "Consider using: ", paste(encas12a_methods, collapse = ", "))
  }

  # ---- Validate DeepHF variant parameter usage ----
  if (score_method != "deephf" && !missing(deephf_var)) {
    message("Note: 'deephf_var' parameter is only used with score_method='deephf'. Ignoring.")
  }

  ## ---- 0. Early‑exit for intragenic mode (single‑exon and two‑exon transcript edge cases) ----------
  n_exons <- length(exon_gr)
  if (n_exons <= 2) {
    warning("Transcript has ", n_exons, " exon(s). Skipping phase‑compatibility checks.")

    # Find sites based on nuclease
    grna_sites <- switch(nuclease,
                         "Cas9" = find_cas9_sites(exon_gr, genome),
                         "Cas12a" = find_cas12a_sites(exon_gr, genome, pam = "TTTV"),
                         "enCas12a" = find_cas12a_sites(exon_gr, genome, pam = "TTTN")
    )

    if (is.null(grna_sites) || length(grna_sites) == 0) {
      warning("No target sites found.")
      return(NULL)
    }

    grna_sites <- score_grnas(grna_sites,
                              method = score_method,
                              tracr = tracr,
                              deephf_var = deephf_var)
    return(grna_sites)
  }

  ## ---- 1. Exon info and pair classification ----
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

  ## ---- 2. Collect and score gRNA sites ----
  if (!is.null(scored_grnas)) {
    # Use pre-scored gRNAs if provided
    grna_sites <- scored_grnas
  } else {
    # Find sites based on nuclease
    grna_sites <- switch(nuclease,
                         "Cas9" = find_cas9_sites(exon_gr, genome),
                         "Cas12a" = find_cas12a_sites(exon_gr, genome, pam = "TTTV"),
                         "enCas12a" = find_cas12a_sites(exon_gr, genome, pam = "TTTN")
    )

    if (is.null(grna_sites) || length(grna_sites) == 0) {
      warning("No sites found within provided exons.")
      return(NULL)
    }

    grna_sites <- score_grnas(grna_sites,
                              method = score_method,
                              tracr = tracr,
                              deephf_var = deephf_var)
  }

  ## ---- 3. Annotate exon rank if missing ----
  if (!"exon_rank" %in% names(mcols(grna_sites))) {
    hits <- GenomicRanges::findOverlaps(grna_sites, exon_gr)
    grna_sites$exon_rank <- NA_integer_
    grna_sites$exon_rank[queryHits(hits)] <- exon_meta$rank[subjectHits(hits)]
  }

  ## ---- 4. Filter to allowed exon ranks ----
  valid <- grna_sites[grna_sites$exon_rank %in% accepted_exons]

  message("Retained ", length(valid),
          " scored gRNAs within phase‑compatible or terminal‑exon‑PTC cases.")

  valid
}
