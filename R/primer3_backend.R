#' Internal function to run Primer3 via reticulate (Batch Mode)
#'
#' Defines and calls a Python function to design primer pairs for a batch of inputs.
#' This offers orders-of-magnitude performance improvements over iterative calls
#' by minimizing R-to-Python context switching.
#'
#' @param request_list List of named lists. Each element must contain:
#'        sequence_template, target_start, target_len, tm_opt, prod_min, prod_max.
#'
#' @return A list of results (each element is a list with keys fwd_seq, rev_seq, prod_size, or NULL).
#' @noRd
run_primer3_batch <- function(request_list) {

  # --- 1. Environment Checks ---
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("reticulate missing")
  if (length(request_list) == 0) return(list())

  if (!check_mutater_env()) {
    if (!activate_mutater_env()) {
      stop("Could not activate mutateR python environment.")
    }
  }

  # --- 2. Define Python Batch Function ---
  is_defined <- tryCatch({
    reticulate::py_eval("'design_primers_batch' in globals()")
  }, error = function(e) FALSE)

  if (!is_defined) {
    reticulate::py_run_string("
import primer3

def design_primers_batch(request_list):
    results = []

    # Define strict global constraints once
    # We copy this dict for every request
    base_global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 65.0,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_GC_CLAMP': 1,         # Weighted preference for 3' GC clamp
        'PRIMER_MAX_SELF_ANY': 8.00,
        'PRIMER_MAX_SELF_END': 3.00,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_NUM_RETURN': 1
    }

    for req in request_list:
        try:
            # Extract request-specific parameters
            seq = req['sequence_template']
            t_start = int(req['target_start'])
            t_len = int(req['target_len'])
            tm_opt = float(req['tm_opt'])
            p_min = int(req['prod_min'])
            p_max = int(req['prod_max'])

            seq_args = {
                'SEQUENCE_ID': 'mutateR_batch',
                'SEQUENCE_TEMPLATE': seq,
                'SEQUENCE_TARGET': [t_start, t_len]
            }

            # Update globals with dynamic Tm/Size settings for this specific request
            current_globals = base_global_args.copy()
            current_globals.update({
                'PRIMER_PRODUCT_SIZE_RANGE': [[p_min, p_max]],
                'PRIMER_OPT_TM': tm_opt,
                'PRIMER_MIN_TM': tm_opt - 5.0,
                'PRIMER_MAX_TM': tm_opt + 5.0
            })

            # Run Primer3
            res = primer3.bindings.design_primers(seq_args, current_globals)

            if res['PRIMER_PAIR_NUM_RETURNED'] > 0:
                results.append({
                    'fwd_seq': res['PRIMER_LEFT_0_SEQUENCE'],
                    'rev_seq': res['PRIMER_RIGHT_0_SEQUENCE'],
                    'prod_size': res['PRIMER_PAIR_0_PRODUCT_SIZE']
                })
            else:
                results.append(None)

        except Exception as e:
            # On error, append None so index alignment is preserved
            # print(f'Batch Error: {e}')
            results.append(None)

    return results
")
  }

  # --- 3. Run Batch Inference ---
  # reticulate automatically converts the R list of lists into a Python list of dicts
  results <- reticulate::py$design_primers_batch(request_list)

  return(results)
}
