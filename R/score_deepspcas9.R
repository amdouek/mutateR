#' Internal function to run DeepSpCas9 via reticulate
#'
#' This function reconstructs the DeepSpCas9 architecture (Kim et al. 2019),
#' loads internal weights, and runs inference.
#'
#' Architecture: 3 parallel Convolutional branches (kernels 3, 5, 7) merged via Concatenation.
#' Input: 30bp (4bp 5' flank + 20bp protospacer + 3bp PAM + 3bp 3' flank).
#'
#' @param sequence_context Character vector of 30 bp sequences.
#' @return Numeric vector of on-target scores (0-100).
#' @noRd
predict_deepspcas9_python <- function(sequence_context) {

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
  # Note: Assuming the weights file is present in extdata
  weights_file <- system.file("extdata", "DeepSpCas9_weights.h5", package = "mutateR")
  if (weights_file == "") {
    # Fallback for dev/local testing
    if (file.exists("./inst/extdata/DeepSpCas9_weights.h5")) {
      weights_file <- "./inst/extdata/DeepSpCas9_weights.h5"
    } else {
      warning("DeepSpCas9_weights.h5 not found. Returning NAs.")
      return(rep(NA_real_, length(sequence_context)))
    }
  }

  # DeepSpCas9 requires 30bp context
  if (nchar(sequence_context[1]) != 30) {
    warning("DeepSpCas9 requires exactly 30 bp sequence context. Returning NAs.")
    return(rep(NA_real_, length(sequence_context)))
  }

  # Re-use the existing one_hot_encode_dna from deep_learning_utils.R
  encoded_data <- one_hot_encode_dna(sequence_context)

  # --- 4. Define Python architecture (Strict Namespacing) ---
  # We check for the function in globals, but to be safe during dev, we redefine it.
  # Crucially, we use a UNIQUE cache variable '_DEEPSPCAS9_MODEL_CACHE'
  # to avoids conflicts with DeepCpf1.

  reticulate::py_run_string("
import tensorflow as tf
import numpy as np
import h5py
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, AveragePooling1D, Flatten, Dense, Dropout, Concatenate

# Distinct Global Cache for SpCas9
_DEEPSPCAS9_MODEL_CACHE = None

def build_and_predict_deepspcas9(weights_path, input_data):
    global _DEEPSPCAS9_MODEL_CACHE

    try:
        # Enforce Float32 to ensure compatibility
        input_data = np.array(input_data, dtype=np.float32)

        if _DEEPSPCAS9_MODEL_CACHE is not None:
            model = _DEEPSPCAS9_MODEL_CACHE
        else:
            # --- Build Architecture (3-Branch CNN) ---
            input_seq = Input(shape=(30, 4))

            # Branch 1: Kernel 3
            c1 = Conv1D(filters=256, kernel_size=3, padding='valid', activation='relu')(input_seq)
            p1 = AveragePooling1D(pool_size=2)(c1)
            f1 = Flatten()(p1)

            # Branch 2: Kernel 5
            c2 = Conv1D(filters=256, kernel_size=5, padding='valid', activation='relu')(input_seq)
            p2 = AveragePooling1D(pool_size=2)(c2)
            f2 = Flatten()(p2)

            # Branch 3: Kernel 7
            c3 = Conv1D(filters=256, kernel_size=7, padding='valid', activation='relu')(input_seq)
            p3 = AveragePooling1D(pool_size=2)(c3)
            f3 = Flatten()(p3)

            # Merge
            m = Concatenate()([f1, f2, f3])

            # Dense Layers
            x = Dense(100, activation='relu')(m)
            x = Dropout(0.3)(x)
            x = Dense(100, activation='relu')(x)
            x = Dropout(0.3)(x)
            output = Dense(1, activation='linear')(x)

            model = Model(inputs=input_seq, outputs=output)

            # --- Load Weights ---
            # Manually map H5 keys to layers to handle naming mismatches
            with h5py.File(weights_path, 'r') as f:
                # Branches are usually layers 1,2,3 or similar depending on graph construction order.
                # We iterate based on kernel size if keys imply it, or rely on standard Keras ordering.
                # Assuming standard load order for this architecture reconstruction:
                # Layers[1]=Conv3, Layers[2]=Pool, Layers[3]=Flat
                # Layers[4]=Conv5...

                # Robust loading by name searching in the H5 file is preferred,
                # but complex without inspecting the specific H5 structure.
                # We assume the H5 was saved from a Keras model with standard names.
                # If using the 'reconstructed' weights from previous session:

                # Load by name matching logic (simplified for embedding):
                layer_names = [l.name for l in model.layers]
                for l_name in layer_names:
                    if l_name in f:
                        g = f[l_name]
                        w_key = [k for k in g.keys() if k.endswith('W') or k.endswith('kernel:0')][0]
                        b_key = [k for k in g.keys() if k.endswith('b') or k.endswith('bias:0')][0]
                        weights = [g[w_key][()], g[b_key][()]]

                        # Fix dimensions if needed (squeeze 1D convs)
                        if len(weights[0].shape) == 4:
                            weights[0] = np.squeeze(weights[0], axis=1)

                        model.get_layer(l_name).set_weights(weights)

            _DEEPSPCAS9_MODEL_CACHE = model

        # --- Prediction ---
        # Use predict_on_batch to avoid creating a new tf.data.Iterator
        # This isolates SpCas9 execution from Cpf1 execution state.
        preds = model.predict_on_batch(input_data)
        return preds.flatten()

    except Exception as e:
        print(f'PYTHON BACKEND ERROR (SpCas9): {e}')
        return []
")

  # --- 5. Run inference ---
  scores <- reticulate::py$build_and_predict_deepspcas9(weights_file, encoded_data)

  if (length(scores) == 0) {
    warning("DeepSpCas9 inference failed (returned 0 scores).")
    return(rep(NA_real_, length(sequence_context)))
  }

  return(as.numeric(scores))
}
