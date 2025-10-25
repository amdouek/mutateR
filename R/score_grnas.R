#' Score individual gRNAs using crisprScore models
#'
#' Handles RuleSet1 (Doench et al., Nat Biotech (2014)) and Azimuth (Doench et al. Nat Biotech (2016) -- in progress) with proper PAM alignment.
#' For models unsupported on Windows, returns NA scores and a warning.
#'
#' @param grna_gr GRanges returned by find_cas9_sites() or find_cas12a_sites().
#' @param method Scoring model name.  Default "ruleset1".
#' @return GRanges with added metadata column `ontarget_score` (numeric)
#' @export
score_grnas <- function(grna_gr,
                        method = c("ruleset1","azimuth","ruleset3",
                                   "deepspcas9","deephf",
                                   "deepcpf1","enpamgb")) {
  
  method <- match.arg(method)
  os_is_windows <- identical(.Platform$OS.type, "windows")
  
  # ---- guard for unsupported models on Windows ----------------------------
  unsupported_windows <- c("ruleset3","deepspcas9","deephf","deepcpf1","enpamgb")
  if (os_is_windows && tolower(method) %in% unsupported_windows) {
    warning("The ", method,
            " model is not supported on Windows. Returning NA scores.")
    mcols(grna_gr)$ontarget_score <- rep(NA_real_, length(grna_gr))
    return(grna_gr)
  }
  
  # ---- basic validation ---------------------------------------------------
  if (!inherits(grna_gr, "GRanges"))
    stop("Input must be a GRanges object.")
  if (!"sequence_context" %in% names(mcols(grna_gr)))
    stop("GRanges must contain 'sequence_context' column.")
  if (!requireNamespace("crisprScore", quietly = TRUE))
    stop("Please install 'crisprScore'.")
  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Please install 'Biostrings'.")
  
  seqs <- as.character(mcols(grna_gr)$sequence_context)
  strand_info <- as.character(strand(grna_gr))
  
  # Reverse‑complement minus‑strand sequences
  is_minus <- strand_info == "-"
  if (any(is_minus)) {
    seqs[is_minus] <- vapply(
      seqs[is_minus],
      \(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))),
      character(1)
    )
  }
  
  # Calculate GC content
  gc_content <- function(s)
    mean(unlist(strsplit(toupper(s), "")) %in% c("G","C"))
  mcols(grna_gr)$gc <- vapply(seqs, gc_content, numeric(1))
  
  message("Computing on‑target scores using ", method, " model...")
  
  # Helper function: Returns numeric vector from any object
  extract_numeric_scores <- function(x) {
    if (is.null(x)) return(rep(NA_real_, length(seqs)))
    if (is.numeric(x)) return(as.numeric(x))
    if (is.data.frame(x) && "score" %in% names(x)) return(as.numeric(x$score))
    if (is.data.frame(x) && "scores" %in% names(x)) return(as.numeric(x$scores))
    if (is.list(x) && "score" %in% names(x)) return(as.numeric(unlist(x$score)))
    if (inherits(x, "DataFrame") && "score" %in% names(x)) return(as.numeric(x$score))
    suppressWarnings(as.numeric(unlist(x)))
  }
  
  scores_raw <- NULL
  
  # ---- SpCas9 30‑bp window models ----------------------------------------
  if (tolower(method) == "ruleset1") {
    not_ngg <- substring(seqs, 26, 27) != "GG"
    if (any(not_ngg, na.rm = TRUE)) {
      message("Re‑orienting ", sum(not_ngg, na.rm = TRUE),
              " guides so PAM = NGG at positions 26–27.")
      seqs[not_ngg] <- vapply(
        seqs[not_ngg],
        \(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))),
        character(1)
      )
    }
    pam_summary <- table(substring(seqs, 25, 27))
    message("PAM triplet distribution (positions 25–27):")
    print(pam_summary)
    
    scores_raw <- crisprScore::getRuleSet1Scores(toupper(seqs))
  }
  
  if (tolower(method) == "azimuth") {
    # enforce PAM = NGG context
    not_ngg <- substring(seqs, 26, 27) != "GG"
    if (any(not_ngg, na.rm = TRUE)) {
      message("Re‑orienting ", sum(not_ngg, na.rm = TRUE),
              " guides so PAM = NGG at positions 26–27 (Azimuth requirement).")
      seqs[not_ngg] <- vapply(
        seqs[not_ngg],
        \(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))),
        character(1)
      )
    }
    pam_summary <- table(substring(seqs, 25, 27))
    message("PAM triplet distribution (positions 25–27):")
    print(pam_summary)
    
    # Skip any remaining non‑NGG guides for model stability
    valid_idx <- substring(seqs, 26, 27) == "GG"
    scores_vec <- rep(NA_real_, length(seqs))
    if (any(valid_idx)) {
      az <- crisprScore::getAzimuthScores(toupper(seqs[valid_idx]))
      scores_vec[valid_idx] <- extract_numeric_scores(az)
    }
    scores_raw <- scores_vec
  }
  
  # ---- DeepHF (Unix only) fallback ---------------------------------------
  if (tolower(method) == "deephf" && !os_is_windows) {
    seq23 <- substring(seqs, 5, 27)
    scores_raw <- crisprScore::getDeepHFScores(toupper(seq23))
  }
  
  # ---- Cas12a models (Unix only) -----------------------------------------
  if (tolower(method) %in% c("deepcpf1","enpamgb") && !os_is_windows) {
    if (tolower(method) == "deepcpf1")
      scores_raw <- crisprScore::getDeepCpf1Scores(toupper(seqs))
    if (tolower(method) == "enpamgb")
      scores_raw <- crisprScore::getEnPAMGBScores(toupper(seqs))
  }
  
  scores_num <- extract_numeric_scores(scores_raw)
  
  mcols(grna_gr)$ontarget_score <- scores_num
  message("Scored ", length(grna_gr), " guides using ", method, ".")
  grna_gr
}