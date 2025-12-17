#' Internal function to run RuleSet3 (Sequence) via reticulate
#'
#' Calls the rs3 Python package (DeWeirdt et al., 2022) for on-target scoring.
#' Uses the Sequence-only model, which predicts activity based on:
#' - Nucleotide identity and position
#' - GC content and thermodynamic properties
#' - tracrRNA scaffold variant ("Chen2013" or "Hsu2013")
#'
#' @param sequence_context Character vector of 30 bp sequences
#'        (4bp 5' flank + 20bp protospacer + 3bp PAM + 3bp 3' flank).
#' @param tracr Character. tracrRNA scaffold: "Chen2013" (default) or "Hsu2013".
#' @param n_jobs Integer. Parallel jobs for featurisation (default 1).
#'
#' @return Numeric vector of RS3 Sequence scores.
#'
#' @references
#' DeWeirdt, P.C., McGee, A.V., Zheng, F., Nwolah, I., Hegde, M. and Doench, J.G., 2022.
#' Accounting for small variations in the tracrRNA sequence improves sgRNA activity
#' predictions for CRISPR screening.
#' Nature Communications, 13(1), p.5255. \doi{10.1038/s41467-022-33024-2}
#'
#' @noRd
predict_rs3_python <- function(sequence_context,
                               tracr = "Chen2013",
                               n_jobs = 1L) {



  # --- 1. Input validation ---

  if (is.null(sequence_context) || length(sequence_context) == 0) {
    return(numeric(0))
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate package is required.")
  }

  # Validate tracr

  valid_tracrs <- c("Hsu2013", "Chen2013")
  if (!tracr %in% valid_tracrs) {
    stop("tracr must be one of: ", paste(valid_tracrs, collapse = ", "))
  }

  # Validate sequence lengths
  seq_lengths <- nchar(sequence_context)
  if (any(seq_lengths != 30)) {
    bad_idx <- which(seq_lengths != 30)
    warning("RS3 requires exactly 30 bp sequence context. ",
            length(bad_idx), " sequences have incorrect length. Returning NAs for those.")
  }


  # --- 2. Environment setup ---

  if (!check_mutater_env()) {
    if (!activate_mutater_env()) {
      stop("Could not activate mutateR Python environment. Please run install_mutater_env().")
    }
  }

  # Check rs3 module availability
  if (!reticulate::py_module_available("rs3")) {
    stop("The 'rs3' Python module is not available.\n",
         "Please reinstall the environment: install_mutater_env(fresh = TRUE)")
  }


  # --- 3. Define Python wrapper function ---

  is_defined <- tryCatch({
    reticulate::py_eval("'predict_rs3_seq' in globals()")
  }, error = function(e) FALSE)

  if (!is_defined) {
    reticulate::py_run_string("
import numpy as np
from rs3.seq import predict_seq

def predict_rs3_seq(sequences, tracr_variant, n_jobs):
    '''
    Wrapper for rs3.seq.predict_seq with error handling.

    Parameters
    ----------
    sequences : list of str
        30-mer context sequences
    tracr_variant : str
        'Hsu2013' or 'Chen2013'
    n_jobs : int
        Parallel jobs for featurization

    Returns
    -------
    numpy.ndarray
        RS3 Sequence scores
    '''
    try:
        # Filter valid sequences (must be exactly 30 bp)
        valid_mask = [len(s) == 30 for s in sequences]
        valid_seqs = [s for s, v in zip(sequences, valid_mask) if v]

        if len(valid_seqs) == 0:
            return np.array([np.nan] * len(sequences))

        # Run RS3 prediction
        scores = predict_seq(
            context_sequences=valid_seqs,
            sequence_tracr=tracr_variant,
            n_jobs=int(n_jobs)
        )

        # Reconstruct full array with NaNs for invalid sequences
        if len(valid_seqs) == len(sequences):
            return np.array(scores)
        else:
            full_scores = np.full(len(sequences), np.nan)
            valid_idx = 0
            for i, v in enumerate(valid_mask):
                if v:
                    full_scores[i] = scores[valid_idx]
                    valid_idx += 1
            return full_scores

    except Exception as e:
        print(f'RS3 PYTHON BACKEND ERROR: {e}')
        return np.array([])
")
  }


  # --- 4. Run inference ---

  scores <- reticulate::py$predict_rs3_seq(
    sequences = as.list(sequence_context),
    tracr_variant = tracr,
    n_jobs = as.integer(n_jobs)
  )

  if (length(scores) == 0) {
    warning("RS3 inference failed (returned 0 scores).")
    return(rep(NA_real_, length(sequence_context)))
  }

  return(as.numeric(scores))
}


#' Validate RS3 installation
#'
#' Tests that the RS3 Python module is correctly installed and functional.
#'
#' @return Logical; TRUE if RS3 is working correctly, FALSE otherwise.
#' @export
validate_rs3 <- function() {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    message("reticulate package not available.")
    return(FALSE)
  }

  if (!check_mutater_env()) {
    if (!activate_mutater_env()) {
      message("Could not activate mutateR environment.")
      return(FALSE)
    }
  }

  if (!reticulate::py_module_available("rs3")) {
    message("RS3 module not found in current Python environment.")
    return(FALSE)

  }

  # Test with a known sequence
  tryCatch({
    test_seq <- "GACGAAAGCGACAACGCGTTCATCCGGGCA"  # From rs3 test data
    scores <- predict_rs3_python(test_seq, tracr = "Chen2013")

    if (length(scores) == 1 && is.numeric(scores) && !is.na(scores)) {
      message("RS3 installation validated successfully.")
      message(sprintf("Test score (Chen2013 tracrRNA): %.4f", scores[1]))
      return(TRUE)
    } else {
      message("RS3 returned unexpected output.")
      return(FALSE)
    }
  }, error = function(e) {
    message(sprintf("RS3 validation failed: %s", e$message))
    return(FALSE)
  })
}
