#' @title On-target scoring for gRNAs
#'
#' @description Performs on-target scoring for gRNAs using selected methods for the relevant nuclease.
#' Currently handles RuleSet1 and Azimuth (Cas9), DeepSpCas9 (Cas9), RuleSet 3 (Cas9),
#' DeepHF (Cas9/eSpCas9/SpCas9-HF1), DeepCpf1 (Cas12a), and enPAM+GB (Cas12a/enCas12a).
#'
#' @param grna_gr GRanges. Returned by find_cas9_sites() or find_cas12a_sites().
#' @param method Character. Scoring model name. Default "ruleset1".
#' @param tracr Character. For Rule Set 3 scoring; one of "Chen2013" (default) or "Hsu2013".
#' @param deephf_var Character. For DeepHF scoring; one of "wt", "wt_u6" (default), "wt_t7", "esp", or "hf".
#'        See \code{\link{recommend_deephf_model}} for guidance on model selection.
#'
#' @section DeepHF Model Selection:
#' The DeepHF method supports multiple enzyme/context variants:
#' \describe{
#'   \item{wt}{Original DeepWt model (base training)}
#'   \item{wt_u6}{Optimised for U6 promoter-driven gRNA expression (default).
#'               Best for: plasmid transfection, lentiviral delivery.}
#'   \item{wt_t7}{Optimised for T7 promoter-driven in vitro transcribed gRNAs.
#'               Best for: RNP delivery with IVT gRNA, synthetic gRNA.}
#'   \item{esp}{Model for eSpCas9(1.1) high-fidelity variant.}
#'   \item{hf}{Model for SpCas9-HF1 high-fidelity variant.}
#' }
#'
#' @return GRanges with added metadata columns `ontarget_score` and `scoring_method`.
#' @export
score_grnas <- function(grna_gr,
                        method = c("ruleset1", "azimuth", "ruleset3",
                                   "deepspcas9", "deephf",
                                   "deepcpf1", "enpamgb"),
                        tracr = "Chen2013",
                        deephf_var = c("wt_u6", "wt_t7", "wt", "esp", "hf")) {

  method <- match.arg(method)
  deephf_var <- match.arg(deephf_var)

  # ---- 1. Basic Validation ----
  if (!inherits(grna_gr, "GRanges")) stop("Input must be a GRanges object.")
  if (!"sequence_context" %in% names(mcols(grna_gr))) stop("GRanges must contain 'sequence_context' column.")

  seqs <- as.character(mcols(grna_gr)$sequence_context)

  # Calculate GC content
  gc_content <- function(s) mean(unlist(strsplit(toupper(s), "")) %in% c("G", "C"))
  mcols(grna_gr)$gc <- vapply(seqs, gc_content, numeric(1))

  message("Computing on-target scores using ", method, " model...")

  # Helper to extract numeric scores
  extract_numeric_scores <- function(x) {
    if (is.null(x)) return(rep(NA_real_, length(seqs)))
    if (is.numeric(x)) return(as.numeric(x))
    if (is.data.frame(x) && "score" %in% names(x)) return(as.numeric(x$score))
    if (is.data.frame(x) && "scores" %in% names(x)) return(as.numeric(x$scores))
    suppressWarnings(as.numeric(unlist(x)))
  }

  scores_raw <- NULL

  # ---- 2. SpCas9 Models (Probability-based) ----
  # Rule Set 1
  if (tolower(method) == "ruleset1") {
    if (!requireNamespace("crisprScore", quietly = TRUE)) stop("crisprScore required.")
    scores_raw <- crisprScore::getRuleSet1Scores(toupper(seqs))
  }

  # Azimuth - not yet tested for functionality - TO DO
  if (tolower(method) == "azimuth") {
    if (!requireNamespace("crisprScore", quietly = TRUE)) stop("crisprScore required.")
    valid_idx <- substring(seqs, 26, 27) == "GG"
    scores_vec <- rep(NA_real_, length(seqs))
    if (any(valid_idx)) {
      az <- crisprScore::getAzimuthScores(toupper(seqs[valid_idx]))
      scores_vec[valid_idx] <- extract_numeric_scores(az)
    }
    scores_raw <- scores_vec
  }

  # Rule Set 3 (via Python backend)
  if (tolower(method) == "ruleset3") {
    message("Using mutateR internal Python backend for RS3 (tracrRNA: ", tracr, ")...")

    if (!check_mutater_env()) {
      warning("Python environment not ready. Running install_mutater_env()...")
      install_mutater_env()
    }

    scores_raw <- tryCatch({
      predict_rs3_python(seqs, tracr = tracr, n_jobs = 1L)
    }, error = function(e) {
      warning("RS3 inference failed: ", e$message)
      return(rep(NA_real_, length(seqs)))
    })
  }

  # ---- 3. DeepSpCas9 (via Python backend) ----
  if (tolower(method) == "deepspcas9") {
    message("Using mutateR internal Python backend for DeepSpCas9...")

    if (!check_mutater_env()) {
      warning("Python environment not ready. Running install_mutater_env()...")
      install_mutater_env()
    }

    scores_raw <- tryCatch({
      predict_deepspcas9_python(seqs)
    }, error = function(e) {
      warning("DeepSpCas9 inference failed: ", e$message)
      return(rep(NA_real_, length(seqs)))
    })
  }

  # ---- 4. DeepHF (via Python backend) ----
  if (tolower(method) == "deephf") {
    message("Using mutateR internal Python backend for DeepHF (variant: ", deephf_var, ")...")

    if (!check_mutater_env()) {
      warning("Python environment not ready. Running install_mutater_env()...")
      install_mutater_env()
    }

    # Check ViennaRNA availability
    if (!reticulate::py_module_available("RNA")) {
      warning("ViennaRNA module not available. DeepHF requires ViennaRNA for RNA structure features.\n",
              "Please reinstall environment: install_mutater_env(fresh = TRUE)")
      scores_raw <- rep(NA_real_, length(seqs))
    } else {
      scores_raw <- tryCatch({
        predict_deephf_python(seqs, deephf_var = deephf_var)
      }, error = function(e) {
        warning("DeepHF inference failed: ", e$message)
        return(rep(NA_real_, length(seqs)))
      })
    }
  }

  # ---- 5. Cas12a Models ----
  # DeepCpf1 (via Python backend)
  if (tolower(method) == "deepcpf1") {
    message("Using mutateR internal Python backend for DeepCpf1...")

    if (!check_mutater_env()) {
      warning("Python environment not ready. Running install_mutater_env()...")
      install_mutater_env()
    }

    scores_raw <- tryCatch({
      predict_deepcpf1_python(seqs)
    }, error = function(e) {
      warning("DeepCpf1 inference failed: ", e$message)
      return(rep(NA_real_, length(seqs)))
    })
  }

  # enPAM+GB (via Python backend - sklearn-based)
  if (tolower(method) == "enpamgb") {
    message("Using mutateR internal Python backend for enPAM+GB...")

    if (!check_mutater_env()) {
      warning("Python environment not ready. Running install_mutater_env()...")
      install_mutater_env()
    }

    scores_raw <- tryCatch({
      predict_enpamgb_python(seqs)
    }, error = function(e) {
      warning("enPAM+GB inference failed: ", e$message)
      return(rep(NA_real_, length(seqs)))
    })
  }

  # ---- 6. Final Formatting ----
  mcols(grna_gr)$ontarget_score <- extract_numeric_scores(scores_raw)
  mcols(grna_gr)$scoring_method <- method

  # Add variant info for DeepHF
  if (tolower(method) == "deephf") {
    mcols(grna_gr)$deephf_var <- deephf_var
  }

  message("Scored ", length(grna_gr), " guides using ", method, ".")

  return(grna_gr)
}
