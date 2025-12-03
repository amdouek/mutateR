#' Internal function to run DeepSpCas9 via reticulate
#'
#' This function reconstructs the DeepSpCas9 architecture (Kim et al. 2019)
#' based on the verified checkpoint structure:
#' - Input: 30bp (4bp 5' flank + 20bp protospacer + 3bp PAM + 3bp 3' flank).
#' - Conv Layer: Inception Module (K3:100, K5:70, K7:40).
#' - Pooling: AveragePooling1D (pool_size=2) on all branches.
#' - Merge: Concatenation (Flattened). Total features = 2790.
#' - Dense: 80 -> 60 -> 1.
#'
#' @param sequence_context Character vector of 30 bp sequences.
#' @return Numeric vector of on-target scores (Linear regression output).
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
  # We look for the converted H5 file now
  weights_file <- system.file("extdata", "DeepSpCas9_weights.h5", package = "mutateR")
  if (weights_file == "") {
    if (file.exists("./inst/extdata/DeepSpCas9_weights.h5")) {
      weights_file <- "./inst/extdata/DeepSpCas9_weights.h5"
    } else {
      warning("DeepSpCas9_weights.h5 not found. Please run conversion script using TensorFlow checkpoint files.")
      return(rep(NA_real_, length(sequence_context)))
    }
  }

  if (nchar(sequence_context[1]) != 30) {
    warning("DeepSpCas9 requires exactly 30 bp sequence context. Returning NAs.")
    return(rep(NA_real_, length(sequence_context)))
  }

  # One-hot encode (uses shared utility)
  encoded_data <- one_hot_encode_dna(sequence_context)

  # --- 4. Define Python Architecture ---
  # Only define if not present in globals
  is_defined <- tryCatch({
    reticulate::py_eval("'build_and_predict_deepspcas9' in globals()")
  }, error = function(e) FALSE)

  if (!is_defined) {
    reticulate::py_run_string("
import tensorflow as tf
import numpy as np
import h5py
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, AveragePooling1D, Flatten, Dense, Dropout, Concatenate

_DEEPSPCAS9_MODEL_CACHE = None

def build_and_predict_deepspcas9(weights_path, input_data):
    global _DEEPSPCAS9_MODEL_CACHE

    try:
        input_data = np.array(input_data, dtype=np.float32)

        if _DEEPSPCAS9_MODEL_CACHE is not None:
            model = _DEEPSPCAS9_MODEL_CACHE
        else:
            # --- Build Architecture ---
            input_seq = Input(shape=(30, 4))

            # Branch 1: Kernel 3 (100 Filters)
            c1 = Conv1D(filters=100, kernel_size=3, padding='valid', activation='relu')(input_seq)
            p1 = AveragePooling1D(pool_size=2)(c1)
            f1 = Flatten()(p1)
            # Output size: (30-3+1)/2 * 100 = 14 * 100 = 1400 features

            # Branch 2: Kernel 5 (70 Filters)
            c2 = Conv1D(filters=70, kernel_size=5, padding='valid', activation='relu')(input_seq)
            p2 = AveragePooling1D(pool_size=2)(c2)
            f2 = Flatten()(p2)
            # Output size: (30-5+1)/2 * 70 = 13 * 70 = 910 features

            # Branch 3: Kernel 7 (40 Filters)
            c3 = Conv1D(filters=40, kernel_size=7, padding='valid', activation='relu')(input_seq)
            p3 = AveragePooling1D(pool_size=2)(c3)
            f3 = Flatten()(p3)
            # Output size: (30-7+1)/2 * 40 = 12 * 40 = 480 features

            # Merge: Order must be 3 -> 5 -> 7
            m = Concatenate()([f1, f2, f3])
            # Total: 1400 + 910 + 480 = 2790 features (Matches checkpoint Dense input)

            # Dense Layers
            x = Dense(80, activation='relu')(m)
            x = Dropout(0.3)(x)
            x = Dense(60, activation='relu')(x)
            x = Dropout(0.3)(x)
            output = Dense(1, activation='linear')(x)

            model = Model(inputs=input_seq, outputs=output)

            # --- Load Weights (From Clean H5) ---
            with h5py.File(weights_path, 'r') as f:
                # Conv 1
                w1 = f['conv1']['weights'][()]
                b1 = f['conv1']['bias'][()]
                model.layers[1].set_weights([w1, b1]) # c1

                # Conv 2
                w2 = f['conv2']['weights'][()]
                b2 = f['conv2']['bias'][()]
                model.layers[2].set_weights([w2, b2]) # c2

                # Conv 3
                w3 = f['conv3']['weights'][()]
                b3 = f['conv3']['bias'][()]
                model.layers[3].set_weights([w3, b3]) # c3

                # Note: layers[4,5,6] are Pooling layers (no weights)
                # layers[7,8,9] are Flatten (no weights)
                # layer[10] is Concatenate

                # We target Dense layers by finding them or by index if stable.
                # Given strict graph construction, indices should be:
                # 11: Dense 1, 13: Dense 2, 15: Output (skipping Dropout layers 12, 14)

                # Robust approach: iterate layers
                dense_layers = [l for l in model.layers if 'dense' in l.name]

                # Dense 1
                wd1 = f['dense1']['weights'][()]
                bd1 = f['dense1']['bias'][()]
                dense_layers[0].set_weights([wd1, bd1])

                # Dense 2
                wd2 = f['dense2']['weights'][()]
                bd2 = f['dense2']['bias'][()]
                dense_layers[1].set_weights([wd2, bd2])

                # Output
                wo = f['output']['weights'][()]
                bo = f['output']['bias'][()]
                dense_layers[2].set_weights([wo, bo])

            _DEEPSPCAS9_MODEL_CACHE = model

        # --- Prediction ---
        preds = model.predict(input_data, verbose=0)
        return preds.flatten()

    except Exception as e:
        print(f'PYTHON BACKEND ERROR (SpCas9): {e}')
        return []
")
  }

  # --- 5. Run inference ---
  scores <- reticulate::py$build_and_predict_deepspcas9(weights_file, encoded_data)

  if (length(scores) == 0) {
    warning("DeepSpCas9 inference failed.")
    return(rep(NA_real_, length(sequence_context)))
  }

  return(as.numeric(scores))
}
