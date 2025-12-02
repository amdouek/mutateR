#' Check exon pair phase compatibility
#'
#' Determines which pairs of exons are phase-compatible, meaning that
#' if the exons between them were deleted, the original reading frame would be preserved.
#'
#' @param exon_info data.frame or GRanges with columns `rank`, `start_phase`, `end_phase`.
#'        Output from get_exon_structures() can be adapted by renaming "phase" to "start_phase".
#' @param include_contiguous Logical. Whether to include contiguous exon pairs
#'        (default FALSE; only non-contiguous pairs are shown).
#'
#' @return A data.frame of compatible exon pairs with columns:
#'   - `exon_5p` (the upstream exon rank)
#'   - `exon_3p` (the downstream exon rank)
#'   - `end_phase_5p`
#'   - `start_phase_3p`
#'   - `compatible` (TRUE/FALSE)
#'
#' @examples
#' \dontrun{
#' exons <- data.frame(rank=1:5,
#'                     start_phase=c(0,2,0,1,1),
#'                     end_phase=c(2,0,1,1,-1))
#' check_exon_phase(exons)
#' }
#' @export

check_exon_phase <- function(exon_info, include_contiguous = FALSE) {
  # ----- Convert GRanges to data.frame if needed -----
  if (inherits(exon_info, "GRanges")) {
    exon_info <- as.data.frame(mcols(exon_info))
    exon_info$rank <- seq_len(nrow(exon_info))
  }

  # ----- Ensure required columns -----
  required_cols <- c("rank", "start_phase", "end_phase")
  if (!all(required_cols %in% colnames(exon_info))) {
    stop("exon_info must contain columns: rank, start_phase, end_phase")
  }

  exon_ranks <- exon_info$rank
  pairs <- t(combn(exon_ranks, 2))
  colnames(pairs) <- c("exon_5p", "exon_3p")
  pairs <- as.data.frame(pairs)

  # ----- Merge upstream end_phase -----
  merged <- merge(pairs, exon_info[, c("rank","end_phase")],
                  by.x="exon_5p", by.y="rank", all.x=TRUE)
  names(merged)[ncol(merged)] <- "end_phase_5p"

  # ----- Merge downstream start_phase -----
  merged <- merge(merged, exon_info[, c("rank","start_phase")],
                  by.x="exon_3p", by.y="rank", all.x=TRUE)
  names(merged)[ncol(merged)] <- "start_phase_3p"

  # ----- Compatibility check -----
  merged$compatible <- merged$end_phase_5p == merged$start_phase_3p

  # ----- Remove contiguous pairs (as default behaviour) -----
  if (!include_contiguous) {
    merged <- merged[(merged[["exon_3p"]] - merged[["exon_5p"]]) > 1, ]
  }

  # ----- Sort nicely (5' - 3') -----
  merged <- merged[order(merged[["exon_5p"]], merged[["exon_3p"]]), ]

  # ----- Reorder columns explicitly to enforce consistent output format -----
  merged <- merged[, c("exon_5p", "exon_3p", "end_phase_5p", "start_phase_3p", "compatible")]

  rownames(merged) <- NULL
  return(merged)
}
