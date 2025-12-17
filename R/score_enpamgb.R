#' Internal function to run enPAM+GB via reticulate
#'
#' This function loads the enPAM+GB scikit-learn Gradient Boosting model
#' (DeWeirdt et al.) and performs feature engineering for Cas12a/enCas12a
#' on-target activity prediction.
#'
#' The model uses the following features:
#' - GC content of protospacer
#' - Position-independent 1-mer and 2-mer frequencies
#' - Position-dependent 1-mer and 2-mer binary indicators
#' - Melting temperature features (context, guide start/mid/end)
#'
#' Input format (34 bp): [4bp 5' flank][4bp PAM][23bp protospacer][3bp 3' flank]
#'
#' @param sequence_context Character vector of 34 bp sequences.
#' @return Numeric vector of on-target scores (linear regression output).
#'
#' @references
#' DeWeirdt, P. C., Sanson, K. R., Sangree, A. K., Hegde, M., Hanna, R. E.,
#' Feeley, M. N., ... & Doench, J. G. (2021). Optimization of AsCas12a for
#' combinatorial genetic screens in human cells. Nature biotechnology, 39(1), 94-104.
#'
#' @noRd
predict_enpamgb_python <- function(sequence_context) {


  # --- 1. Input validation ---
  if (is.null(sequence_context) || length(sequence_context) == 0) {
    return(numeric(0))
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate package is required.")
  }

  # Validate sequence lengths

  seq_lengths <- nchar(sequence_context)
  if (any(seq_lengths != 34)) {
    bad_idx <- which(seq_lengths != 34)
    warning("enPAM+GB requires exactly 34 bp sequence context. ",
            length(bad_idx), " sequences have incorrect length. Returning NAs for those.")
  }


  # --- 2. Environment setup ---
  if (!check_mutater_env()) {
    if (!activate_mutater_env()) {
      stop("Could not activate mutateR Python environment. Please run install_mutater_env().")
    }
  }

  # Check required modules
  if (!reticulate::py_module_available("sklearn")) {
    stop("scikit-learn is not available in the Python environment.\n",
         "Please reinstall: install_mutater_env(fresh = TRUE)")
  }

  if (!reticulate::py_module_available("Bio")) {
    stop("Biopython is not available in the Python environment.\n",
         "Please reinstall: install_mutater_env(fresh = TRUE)")
  }


  # --- 3. Locate model weights ---
  weights_file <- system.file("extdata", "enPAM_GB.joblib", package = "mutateR")

  if (weights_file == "") {
    # Fallback for dev environments
    if (file.exists("./inst/extdata/enPAM_GB.joblib")) {
      weights_file <- "./inst/extdata/enPAM_GB.joblib"
    } else {
      stop("Expected enPAM+GB weights (enPAM_GB.joblib) not found in package extdata directory.")
    }
  }


  # --- 4. Define Python backend ---
  is_defined <- tryCatch({
    reticulate::py_eval("'predict_enpamgb' in globals()")
  }, error = function(e) FALSE)

  if (!is_defined) {
    reticulate::py_run_string("
import numpy as np
import pandas as pd
import joblib
import sys
import types
import re
from Bio.SeqUtils import MeltingTemp

# ============================================================================
# SKLEARN LEGACY PATCHES
# Required for loading older .joblib files with newer sklearn versions
# ============================================================================

def patch_sklearn_gradient_boosting():
    '''Recreate old sklearn.ensemble.gradient_boosting module path'''
    from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor

    mod_name = 'sklearn.ensemble.gradient_boosting'
    module = types.ModuleType(mod_name)
    module.GradientBoostingClassifier = GradientBoostingClassifier
    module.GradientBoostingRegressor = GradientBoostingRegressor
    sys.modules[mod_name] = module


def patch_sklearn_tree_tree():
    '''Recreate old sklearn.tree.tree module path'''
    from sklearn.tree import (
        DecisionTreeClassifier,
        DecisionTreeRegressor,
        ExtraTreeClassifier,
        ExtraTreeRegressor,
    )

    mod_name = 'sklearn.tree.tree'
    module = types.ModuleType(mod_name)
    module.DecisionTreeClassifier = DecisionTreeClassifier
    module.DecisionTreeRegressor = DecisionTreeRegressor
    module.ExtraTreeClassifier = ExtraTreeClassifier
    module.ExtraTreeRegressor = ExtraTreeRegressor
    sys.modules[mod_name] = module


LEGACY_PATCHES = {
    'sklearn.ensemble.gradient_boosting': patch_sklearn_gradient_boosting,
    'sklearn.tree.tree': patch_sklearn_tree_tree,
}


def load_with_autopatch(path):
    '''Load joblib file with automatic legacy module patching'''
    tried_modules = set()

    while True:
        try:
            return joblib.load(path)
        except ModuleNotFoundError as e:
            msg = str(e)
            m = re.search(r\"No module named '([^']+)'\", msg)
            if not m:
                raise

            missing_module = m.group(1)

            if missing_module in tried_modules:
                raise

            tried_modules.add(missing_module)
            patch_fn = LEGACY_PATCHES.get(missing_module)

            if patch_fn is None:
                raise

            patch_fn()


# ============================================================================
# FEATURE ENGINEERING FUNCTIONS
# Ported from sgrna_modeler/features.py
# ============================================================================

def get_guide_sequence(context, guide_start, guide_length):
    '''Extract guide sequence from context (1-indexed start)'''
    return context[(guide_start - 1):(guide_start + guide_length - 1)]


def get_context_order(k):
    '''Generate position labels for k-mer context'''
    return [str(x + 1) for x in range(k)]


def get_frac_g_or_c(curr_dict, guide_sequence):
    '''Calculate GC content of guide'''
    g_count = guide_sequence.count('G')
    c_count = guide_sequence.count('C')
    gc_frac = (g_count + c_count) / len(guide_sequence)
    curr_dict['GC content'] = gc_frac
    return curr_dict


def get_one_nt_counts(curr_dict, guide, nts):
    '''Position-independent 1-mer frequencies'''
    for nt in nts:
        nt_frac = guide.count(nt) / len(guide)
        curr_dict[nt] = nt_frac
    return curr_dict


def get_two_nt_counts(curr_dict, guide, nts):
    '''Position-independent 2-mer frequencies'''
    for nt1 in nts:
        for nt2 in nts:
            two_mer = nt1 + nt2
            nts_counts = guide.count(two_mer)
            nts_frac = nts_counts / (len(guide) - 1)
            curr_dict[nt1 + nt2] = nts_frac
    return curr_dict


def get_one_nt_pos(curr_dict, context_sequence, nts, context_order):
    '''Position-dependent 1-mer binary indicators'''
    for i in range(len(context_order)):
        curr_nt = context_sequence[i]
        for nt in nts:
            key = context_order[i] + nt
            curr_dict[key] = 1 if curr_nt == nt else 0
    return curr_dict


def get_two_nt_pos(curr_dict, context_sequence, nts, context_order):
    '''Position-dependent 2-mer binary indicators'''
    for i in range(len(context_order) - 1):
        curr_nts = context_sequence[i:i+2]
        for nt1 in nts:
            for nt2 in nts:
                match_nts = nt1 + nt2
                key = context_order[i] + match_nts
                curr_dict[key] = 1 if curr_nts == match_nts else 0
    return curr_dict


def get_thermo(curr_dict, guide_sequence, context_sequence):
    '''Calculate melting temperature features'''
    curr_dict['Tm, context'] = MeltingTemp.Tm_NN(context_sequence)
    third = len(guide_sequence) // 3
    curr_dict['Tm, start'] = MeltingTemp.Tm_NN(guide_sequence[0:third])
    curr_dict['Tm, mid'] = MeltingTemp.Tm_NN(guide_sequence[third:2*third])
    curr_dict['Tm, end'] = MeltingTemp.Tm_NN(guide_sequence[2*third:])
    return curr_dict


def featurize_guides(kmers, guide_start=9, guide_length=23):
    '''
    Featurize guide sequences for enPAM+GB model.

    Features computed:
    - GC content
    - Position-independent 1-mer (4 features)
    - Position-independent 2-mer (16 features)
    - Position-dependent 1-mer (34 * 4 = 136 features)
    - Position-dependent 2-mer (33 * 16 = 528 features)
    - Melting temperature (4 features)

    Total: 689 features

    Parameters
    ----------
    kmers : list of str
        34 bp context sequences
    guide_start : int
        1-indexed start position of guide (9 for Cas12a)
    guide_length : int
        Length of guide sequence (23 for Cas12a)

    Returns
    -------
    pd.DataFrame
        Feature matrix
    '''
    k = len(kmers[0])
    context_order = get_context_order(k)
    nts = ['A', 'C', 'T', 'G']
    feature_dict_list = []

    for i in range(len(kmers)):
        curr_dict = {}
        context = kmers[i]
        guide_sequence = get_guide_sequence(context, guide_start, guide_length)

        # GC content
        curr_dict = get_frac_g_or_c(curr_dict, guide_sequence)

        # Position-independent features
        curr_dict = get_one_nt_counts(curr_dict, guide_sequence, nts)
        curr_dict = get_two_nt_counts(curr_dict, guide_sequence, nts)

        # Position-dependent features
        curr_dict = get_one_nt_pos(curr_dict, context, nts, context_order)
        curr_dict = get_two_nt_pos(curr_dict, context, nts, context_order)

        # Melting temperature features
        curr_dict = get_thermo(curr_dict, guide_sequence, context)

        feature_dict_list.append(curr_dict)

    return pd.DataFrame(feature_dict_list)


# ============================================================================
# MODEL CACHE AND PREDICTION
# ============================================================================

_ENPAMGB_MODEL_CACHE = None


def predict_enpamgb(model_path, sequences):
    '''
    Load enPAM+GB model and predict on-target scores.

    Parameters
    ----------
    model_path : str
        Path to enPAM_GB.joblib file
    sequences : list of str
        34 bp context sequences

    Returns
    -------
    np.ndarray
        Predicted on-target scores
    '''
    global _ENPAMGB_MODEL_CACHE

    try:
        # Validate sequences
        valid_mask = [len(s) == 34 for s in sequences]
        valid_seqs = [s.upper() for s, v in zip(sequences, valid_mask) if v]

        if len(valid_seqs) == 0:
            return np.array([np.nan] * len(sequences))

        # Load model (with caching)
        if _ENPAMGB_MODEL_CACHE is None:
            _ENPAMGB_MODEL_CACHE = load_with_autopatch(model_path)

        model = _ENPAMGB_MODEL_CACHE

        # Featurize sequences
        # Cas12a parameters: guide_start=9, guide_length=23
        feature_matrix = featurize_guides(valid_seqs, guide_start=9, guide_length=23)

        # Predict
        predictions = model.predict(feature_matrix)

        # Reconstruct full array with NaNs for invalid sequences
        if len(valid_seqs) == len(sequences):
            return np.array(predictions)
        else:
            full_scores = np.full(len(sequences), np.nan)
            valid_idx = 0
            for i, v in enumerate(valid_mask):
                if v:
                    full_scores[i] = predictions[valid_idx]
                    valid_idx += 1
            return full_scores

    except Exception as e:
        print(f'enPAM+GB PYTHON BACKEND ERROR: {e}')
        return np.array([])
")
  }


# --- 5. Run inference ---
scores <- reticulate::py$predict_enpamgb(
  model_path = weights_file,
  sequences = as.list(sequence_context)
)

if (length(scores) == 0) {
  warning("enPAM+GB inference failed (returned 0 scores).")
  return(rep(NA_real_, length(sequence_context)))
}

return(as.numeric(scores))
}


#' Validate enPAM+GB installation
#'
#' Tests that the enPAM+GB model and dependencies are correctly installed.
#'
#' @return Logical; TRUE if enPAM+GB is working correctly, FALSE otherwise.
#' @export
validate_enpamgb <- function() {

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

  # Check dependencies
  if (!reticulate::py_module_available("sklearn")) {
    message("scikit-learn not found in current Python environment.")
    return(FALSE)
  }

  if (!reticulate::py_module_available("Bio")) {
    message("Biopython not found in current Python environment.")
    return(FALSE)
  }

  # Test with known sequences from sgrna_modeler test data
  tryCatch({
    test_seqs <- c(
      "TGGTTTTAAAACAGAATATACAGTCTAAAAAACC",
      "CATGTTTTTTTGGGAACCAATCGATAATCACATT"
    )

    scores <- predict_enpamgb_python(test_seqs)

    if (length(scores) == 2 && is.numeric(scores) && !any(is.na(scores))) {
      message("enPAM+GB installation validated successfully.")
      message(sprintf("Test scores: %.4f, %.4f", scores[1], scores[2]))
      return(TRUE)
    } else {
      message("enPAM+GB returned unexpected output.")
      return(FALSE)
    }
  }, error = function(e) {
    message(sprintf("enPAM+GB validation failed: %s", e$message))
    return(FALSE)
  })
}
