#' Assess if an deletion of exon(s) would cause a frameshift/PTC
#'
#' @param exon_info data.frame with `rank` and `exon_cds_length`
#' @param exon_5p Integer. Upstream exon rank.
#' @param exon_3p Integer. Downstream exon rank.
#'
#' @return List with frameshift (logical), ptc_flag (logical),
#'         terminal_exon_case (logical), and deleted_length (integer).
check_frameshift_ptc <- function(exon_info, exon_5p, exon_3p) {
  # ----- Ensure required columns -----
  required <- c("rank", "exon_cds_length")
  if (!all(required %in% names(exon_info))) {
    stop("exon_info must contain: ", paste(required, collapse = ", "))
  }

  # ----- Specify which exons are being deleted -----
  # The exons BETWEEN the flanking pair are removed
  exons_to_delete <- exon_info$rank[(exon_info$rank > exon_5p) & (exon_info$rank < exon_3p)]

  # ----- Calculate total coding bases removed -----
  deleted_length <- sum(exon_info$exon_cds_length[exon_info$rank %in% exons_to_delete],
                        na.rm = TRUE)

  # ----- Determine if exon removal causes frameshift -----
  frameshift <- (deleted_length %% 3) != 0

  # ----- Specify terminal exon as highest rank -----
  terminal_rank <- max(exon_info$rank, na.rm = TRUE)

  ptc_flag <- FALSE
  terminal_exon_case <- FALSE

  if (frameshift) {
    if (exon_3p == terminal_rank) {
      terminal_exon_case <- TRUE
    } else {
      ptc_flag <- TRUE
    }
  }

  return(list(frameshift = frameshift,
              ptc_flag = ptc_flag,
              terminal_exon_case = terminal_exon_case,
              deleted_length = as.integer(deleted_length)))
}
