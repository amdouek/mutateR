#' Score individual gRNAs using crisprScore models
#'
#' Handles RuleSet1 (Doench et al., Nat Biotech (2014)) and Azimuth (Doench et al. Nat Biotech (2016) -- in progress) for Cas9 and DeepCpf1 for Cas12.
#' Automatically routes DeepCpf1 requests to the internal Python backend
#' on Windows to bypass crisprScore OS limitations.
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

  # ---- 1. Unsupported models on Windows (excluding DeepCpf1) ----
  unsupported_windows <- c("ruleset3","deepspcas9","deephf","enpamgb")
  if (os_is_windows && tolower(method) %in% unsupported_windows) {
    warning("The ", method, " model is not supported on Windows. Returning NA scores.")
    mcols(grna_gr)$ontarget_score <- rep(NA_real_, length(grna_gr))
    return(grna_gr)
  }

  # ---- 2. Basic Validation ----------------------------------------
  if (!inherits(grna_gr, "GRanges")) stop("Input must be a GRanges object.")
  if (!"sequence_context" %in% names(mcols(grna_gr))) stop("GRanges must contain 'sequence_context' column.")

  seqs <- as.character(mcols(grna_gr)$sequence_context)

  # ---- RESTORED: PAM Summary Table ----
  if ("pam_sequence" %in% names(mcols(grna_gr))) {
    # If explicit PAM column exists (Cas12a / Updated Cas9)
    pam_summary <- table(mcols(grna_gr)$pam_sequence)
    message("PAM distribution:")
    print(pam_summary)
  } else if (nchar(seqs[1]) >= 27) {
    # Fallback for Cas9 RuleSet1/Azimuth (PAM usually at 25-27 of 30bp context)
    pam_summary <- table(substring(seqs, 25, 27))
    message("PAM distribution (derived from context):")
    print(pam_summary)
  }

  # Calculate GC content
  gc_content <- function(s) mean(unlist(strsplit(toupper(s), "")) %in% c("G","C"))
  mcols(grna_gr)$gc <- vapply(seqs, gc_content, numeric(1))

  message("Computing on‑target scores using ", method, " model...")

  # Helper to extract numeric vector
  extract_numeric_scores <- function(x) {
    if (is.null(x)) return(rep(NA_real_, length(seqs)))
    if (is.numeric(x)) return(as.numeric(x))
    if (is.data.frame(x) && "score" %in% names(x)) return(as.numeric(x$score))
    if (is.data.frame(x) && "scores" %in% names(x)) return(as.numeric(x$scores))
    if (inherits(x, "DataFrame") && "score" %in% names(x)) return(as.numeric(x$score))
    suppressWarnings(as.numeric(unlist(x)))
  }

  scores_raw <- NULL

  # ---- 3. SpCas9 Models -------------------------------------------
  if (tolower(method) == "ruleset1") {
    if (!requireNamespace("crisprScore", quietly = TRUE)) stop("crisprScore required.")
    scores_raw <- crisprScore::getRuleSet1Scores(toupper(seqs))
  }

  if (tolower(method) == "azimuth") {
    if (!requireNamespace("crisprScore", quietly = TRUE)) stop("crisprScore required.")
    # Check for NGG at pos 26-27
    valid_idx <- substring(seqs, 26, 27) == "GG"
    scores_vec <- rep(NA_real_, length(seqs))
    if (any(valid_idx)) {
      az <- crisprScore::getAzimuthScores(toupper(seqs[valid_idx]))
      scores_vec[valid_idx] <- extract_numeric_scores(az)
    }
    scores_raw <- scores_vec
  }

  if (tolower(method) == "deephf" && !os_is_windows) {
    if (!requireNamespace("crisprScore", quietly = TRUE)) stop("crisprScore required.")
    seq23 <- substring(seqs, 5, 27)
    scores_raw <- crisprScore::getDeepHFScores(toupper(seq23))
  }

  # ---- 4. Cas12a Models (DeepCpf1) --------------------------------
  if (tolower(method) == "deepcpf1") {

    # Check context length (Seq-DeepCpf1 needs 34bp)
    if (any(nchar(seqs) != 34)) {
      warning("DeepCpf1 requires 34bp context. Some guides may fail scoring.")
    }

    use_internal_backend <- os_is_windows

    if (!use_internal_backend) {
      # Try crisprScore first on non-Windows
      if (requireNamespace("crisprScore", quietly = TRUE)) {
        tryCatch({
          scores_raw <- crisprScore::getDeepCpf1Scores(toupper(seqs))
        }, error = function(e) {
          message("crisprScore DeepCpf1 failed (", e$message, "). Falling back to internal python.")
          use_internal_backend <<- TRUE
        })
      } else {
        use_internal_backend <- TRUE
      }
    }

    if (use_internal_backend) {
      message("Using mutateR internal Python backend for DeepCpf1...")

      # Ensure environment is ready
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
  }

  # ---- 5. Final Formatting ----------------------------------------
  scores_num <- extract_numeric_scores(scores_raw)
  mcols(grna_gr)$ontarget_score <- scores_num
  message("Scored ", length(grna_gr), " guides using ", method, ".")

  return(grna_gr)
}
