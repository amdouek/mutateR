#' Fetch DeepCpf1 model weights
#'
#' Retrieves the path to the embedded DeepCpf1 weights file.
#' Note: The weights used is for the Seq-DeepCpf1 model (solely sequence-dependent). The better-performing DeepCpf1 model integrates chromatin accessibility data, but this is cell-line specific rather than useful across organisms, and so is not implemented here.
#'
#' @return Path to the weights file (.h5).
#' @export
fetch_deepcpf1_weights <- function() {
  weights_path <- system.file("extdata", "Seq_deepCpf1_weights.h5", package = "mutateR")

  if (weights_path == "") {
    # Fallback for dev environments if needed
    if (file.exists("./inst/extdata/Seq_deepCpf1_weights.h5")) {
      weights_path <- "./inst/extdata/Seq_deepCpf1_weights.h5"
    } else {
      stop("Seq_deepCpf1_weights.h5 not found in package directory.")
    }
  }
  return(weights_path)
}

#' One-hot encode DNA sequences for deep learning
#'
#' Converts a set of DNA sequences into a 3D numeric array
#' (Samples x Length x 4) suitable for Keras/TensorFlow.
#' Optimised for performance.
#'
#' @param sequences Character vector of DNA sequences.
#'
#' @return An array of shape (N, Length, 4).
#' @noRd
one_hot_encode_dna <- function(sequences) {

  if (length(sequences) == 0) return(array(0, dim = c(0, 0, 4)))

  # Ensure uppercase
  sequences <- toupper(sequences)

  # Check lengths
  n_seqs <- length(sequences)
  seq_len <- nchar(sequences[1])

  # Fast matrix construction
  # unlist(strsplit) is much faster than do.call(rbind)
  char_vec <- unlist(strsplit(sequences, "", fixed = TRUE), use.names = FALSE)

  if (length(char_vec) != n_seqs * seq_len) {
    stop("All sequences must be the same length for one-hot encoding.")
  }

  char_mat <- matrix(char_vec, nrow = n_seqs, ncol = seq_len, byrow = TRUE)

  # Initialise 3D array: (N, Length, 4)
  # Channels: A, C, G, T
  x_encoded <- array(0, dim = c(n_seqs, seq_len, 4))

  # Vectorised assignment
  x_encoded[,,1] <- (char_mat == "A") * 1
  x_encoded[,,2] <- (char_mat == "C") * 1
  x_encoded[,,3] <- (char_mat == "G") * 1
  x_encoded[,,4] <- (char_mat == "T") * 1

  return(x_encoded)
}
