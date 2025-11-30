#' Internal function to run DeepCpf1 via reticulate
#'
#' This function defines the Seq-DeepCpf1 architecture (Kim et al. 2018),
#' manually loads the legacy weights from the .h5 file to resolve shape mismatches,
#' and runs inference on the one-hot encoded DNA.
#'
#' @param sequence_context Character vector of 34 bp sequences.
#' @return Numeric vector of on-target scores.
#' @noRd
predict_deepcpf1_python <- function(sequence_context) {

  # --- 1. Input tidying ---
  if (is.null(sequence_context) || length(sequence_context) == 0) {
    return(numeric(0))
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) stop("reticulate missing")

  # --- 2. Environment setup ---
  if (!check_mutater_env()) {
    if (!activate_mutater_env()) {
      stop("Could not activate mutateR python environment. Please run install_mutater_env().")
    }
  }

  # --- 3. Prepare data & weights ---
  weights_file <- fetch_deepcpf1_weights()

  if (nchar(sequence_context[1]) != 34) {
    warning("Seq-DeepCpf1 requires exactly 34 bp sequence context. Returning NAs.")
    return(rep(NA_real_, length(sequence_context)))
  }

  encoded_data <- one_hot_encode_dna(sequence_context)

  # --- 4. Define Python architecture ---
  is_defined <- tryCatch({
    reticulate::py_eval("'build_and_predict_deepcpf1' in globals()")
  }, error = function(e) FALSE)

  # --- Updated Python string with model caching ---
  if (!is_defined) {
    reticulate::py_run_string("
import tensorflow as tf
import numpy as np
import h5py
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, AveragePooling1D, Flatten, Dense, Dropout

# Global cache variable to store the model
_DEEPCPF1_MODEL_CACHE = None

def build_and_predict_deepcpf1(weights_path, input_data):
    global _DEEPCPF1_MODEL_CACHE

    try:
        # Check if we already built the model
        if _DEEPCPF1_MODEL_CACHE is not None:
            model = _DEEPCPF1_MODEL_CACHE
        else:
            # --- Build the Architecture ---
            input_seq = Input(shape=(34, 4))
            x = Conv1D(80, 5, padding='valid', activation='relu')(input_seq)
            x = AveragePooling1D(2)(x)
            x = Flatten()(x)
            x = Dense(80, activation='relu')(x)
            x = Dropout(0.3)(x) # Note: Dropout is ignored during inference
            x = Dense(40, activation='relu')(x)
            x = Dropout(0.3)(x)
            x = Dense(40, activation='relu')(x)
            x = Dropout(0.3)(x)
            output = Dense(1, activation='linear')(x)

            model = Model(inputs=input_seq, outputs=output)

            # --- Load Weights ---
            with h5py.File(weights_path, 'r') as f:
                # Conv Weights
                conv_keys = [k for k in f.keys() if 'convolution1d' in k]
                c_key = conv_keys[0]
                conv_W = f[c_key][f'{c_key}_W'][()]
                conv_b = f[c_key][f'{c_key}_b'][()]
                if len(conv_W.shape) == 4: conv_W = np.squeeze(conv_W, axis=1)
                model.layers[1].set_weights([conv_W, conv_b])

                # Dense Weights
                dense_keys = sorted([k for k in f.keys() if 'dense' in k])
                layer_map = [4, 6, 8, 10]
                for i, d_key in enumerate(dense_keys):
                    if i >= len(layer_map): break
                    target_idx = layer_map[i]
                    dW = f[d_key][f'{d_key}_W'][()]
                    db = f[d_key][f'{d_key}_b'][()]
                    model.layers[target_idx].set_weights([dW, db])

            # Save to cache
            _DEEPCPF1_MODEL_CACHE = model

        # --- Prediction ---
        preds = model.predict(input_data, verbose=0)
        return preds.flatten()

    except Exception as e:
        print(f'PYTHON BACKEND ERROR: {e}')
        return []
")
  }

  # --- 5. Run inference ---
  scores <- reticulate::py$build_and_predict_deepcpf1(weights_file, encoded_data)

  if (length(scores) == 0) {
    warning("DeepCpf1 inference failed (returned 0 scores).")
    return(rep(NA_real_, length(sequence_context)))
  }

  return(as.numeric(scores))
}
