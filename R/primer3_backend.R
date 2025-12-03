#' Internal function to run Primer3 via reticulate
#'
#' Defines and calls a Python function to design primer pairs flanking a specific target region.
#' Implements strict filtering:
#' - Homopolymers (Max 4)
#' - GC Content (40-65%)
#' - 3' GC Clamp (Weighted preference)
#' - Thermodynamic secondary structure checks
#'
#' @param sequence_template Character. The DNA sequence context (WT).
#' @param target_start Integer. 0-based start index of the deletion target within the template.
#' @param target_len Integer. Length of the deletion target.
#' @param tm_opt Numeric. Optimal Tm (default 60.0).
#' @param prod_min Integer. Min product size.
#' @param prod_max Integer. Max product size.
#'
#' @return A list with keys: fwd_seq, rev_seq, product_size (or NULL if failed).
#' @noRd
run_primer3_python <- function(sequence_template,
                               target_start,
                               target_len,
                               tm_opt = 60.0,
                               prod_min = 100,
                               prod_max = 1000) {

  # --- 1. Environment Checks ---
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("reticulate missing")

  if (!check_mutater_env()) {
    if (!activate_mutater_env()) {
      stop("Could not activate mutateR python environment.")
    }
  }

  # --- 2. Define Python Function (if none exists) ---
  is_defined <- tryCatch({
    reticulate::py_eval("'design_primers_python' in globals()")
  }, error = function(e) FALSE)

  if (!is_defined) {
    reticulate::py_run_string("
import primer3

def design_primers_python(seq, t_start, t_len, tm_opt, p_min, p_max):
    try:
        # Define the target region to flank (start, length)
        # Note: primer3 uses 0-based indexing for sequence position
        seq_args = {
            'SEQUENCE_ID': 'mutateR_design',
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_TARGET': [int(t_start), int(t_len)]
        }

        global_args = {
            # --- Size Constraints ---
            'PRIMER_PRODUCT_SIZE_RANGE': [[int(p_min), int(p_max)]],
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,

            # --- Tm Constraints ---
            'PRIMER_OPT_TM': float(tm_opt),
            'PRIMER_MIN_TM': float(tm_opt) - 5.0,
            'PRIMER_MAX_TM': float(tm_opt) + 5.0,

            # --- Composition Constraints ---
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 65.0,
            'PRIMER_MAX_POLY_X': 4,       # Reject runs >4 (e.g. AAAAA)
            'PRIMER_GC_CLAMP': 1,         # 1 = Require 3' GC clamp if possible (penalty otherwise)

            # --- Structure Constraints (Thermodynamic) ---
            # Max Self Alignment (Dimer): < 8.00 score usually safe
            'PRIMER_MAX_SELF_ANY': 8.00,
            # Max 3' End Self Alignment (Primer-Dimer): Strict
            'PRIMER_MAX_SELF_END': 3.00,

            # --- Misc ---
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_NUM_RETURN': 1
        }

        # Updated to snake_case to avoid deprecation warning
        res = primer3.bindings.design_primers(seq_args, global_args)

        # Check if a pair was returned
        if res['PRIMER_PAIR_NUM_RETURNED'] > 0:
            return {
                'fwd_seq': res['PRIMER_LEFT_0_SEQUENCE'],
                'rev_seq': res['PRIMER_RIGHT_0_SEQUENCE'],
                'prod_size': res['PRIMER_PAIR_0_PRODUCT_SIZE']
            }
        return None

    except Exception as e:
        # print(f'Primer3 Error: {e}')
        return None
")
  }

  # --- 3. Run Inference ---
  # Ensure integer types for Python
  res <- reticulate::py$design_primers_python(
    sequence_template,
    as.integer(target_start),
    as.integer(target_len),
    as.numeric(tm_opt),
    as.integer(prod_min),
    as.integer(prod_max)
  )

  return(res)
}
