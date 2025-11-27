#' Interactive gRNA design plot
#'
#' Interactive plotting mode for heatmap via plotly.
#'
#' @export
#' Interactive gRNA design plot
#'
#' Interactive plotting mode for heatmap via plotly.
#' Handles empty pairs_df gracefully by plotting just the phase compatibility.
#'
#' @export
plot_grna_interactive <- function(exon_gr,
                                  pairs_df = NULL,
                                  transcript_id,
                                  gene_symbol = NULL,
                                  species) {

  if (!requireNamespace("plotly", quietly = TRUE)) stop("Package 'plotly' required.")
  library(dplyr)
  library(plotly)

  ## ---- 1. Exon metadata ----
  ex_meta <- as.data.frame(mcols(exon_gr))
  n_exons <- length(exon_gr)
  ex_meta$rank <- seq_len(n_exons)
  ex_meta <- ex_meta %>%
    mutate(is_UTR = (start_phase == -1 & end_phase == -1))

  ## ---- 2. Build complete exon pair grid ----
  all_pairs <- expand.grid(exon_5p = ex_meta$rank,
                           exon_3p = ex_meta$rank,
                           KEEP.OUT.ATTRS = FALSE)

  comp_df <- check_exon_phase(ex_meta, include_contiguous = TRUE)

  merged <- left_join(all_pairs, comp_df, by = c("exon_5p","exon_3p"))

  # Define Categories
  merged <- merged %>%
    mutate(Category = case_when(
      exon_5p == exon_3p         ~ "Self",
      exon_3p - exon_5p == 1     ~ "Contiguous",
      is.na(compatible)          ~ "Incompatible",
      compatible                 ~ "Compatible",
      TRUE                       ~ "Incompatible"
    ))

  utr_exons <- ex_meta %>% filter(is_UTR) %>% pull(rank)
  merged$Category[merged$exon_5p %in% utr_exons |
                    merged$exon_3p %in% utr_exons] <- "UTR"

  ## ---- Color map ----
  cat_levels <- c("Self","Contiguous","Compatible","Incompatible","UTR")
  cat_cols   <- c("grey70","palegreen3","goldenrod2","steelblue4","firebrick3")

  merged$CatNum <- match(merged$Category, cat_levels) - 1

  ## ---- 3. gRNA summary for hover ----
  has_pairs <- !is.null(pairs_df) && nrow(pairs_df) > 0

  if (has_pairs) {
    pair_summary <- pairs_df %>%
      group_by(upstream_pair, downstream_pair) %>%
      summarise(
        cut_5p   = dplyr::first(exon_5p),
        cut_3p   = dplyr::first(exon_3p),
        top5p    = protospacer_sequence_5p[which.max(ontarget_score_5p)],
        top3p    = protospacer_sequence_3p[which.max(ontarget_score_3p)],
        score5p  = max(ontarget_score_5p, na.rm = TRUE),
        score3p  = max(ontarget_score_3p, na.rm = TRUE),
        n_pairs  = n(),
        .groups  = "drop"
      ) %>%
      rename(exon_5p = upstream_pair,
             exon_3p = downstream_pair)

    merged <- left_join(merged, pair_summary, by = c("exon_5p", "exon_3p"))
  }

  ## ---- 4. Construct clean hover labels ----
  if (!"n_pairs" %in% names(merged)) merged$n_pairs <- NA
  if (!"score5p" %in% names(merged)) merged$score5p <- NA
  if (!"score3p" %in% names(merged)) merged$score3p <- NA
  if (!"cut_5p" %in% names(merged))  merged$cut_5p <- NA
  if (!"cut_3p" %in% names(merged))  merged$cut_3p <- NA

  merged <- merged %>%
    mutate(hover = paste0(
      "<b>E", exon_5p, " – E", exon_3p, "</b>",
      "<br>Status: ", Category,
      ifelse(!is.na(n_pairs), paste0("<br>Pairs found: ", n_pairs),
             ifelse(Category == "Compatible", "<br>No valid gRNA pairs", "")),
      ifelse(!is.na(score5p),
             paste0("<br>Top 5' (E", cut_5p, "): ",
                    formatC(score5p, digits = 3, format = 'f')), ""),
      ifelse(!is.na(score3p),
             paste0("<br>Top 3' (E", cut_3p, "): ",
                    formatC(score3p, digits = 3, format = 'f')), "")
    ))

  ## ---- 5. Construct Matrices (Standard Order) ----
  merged <- merged[order(merged$exon_3p, merged$exon_5p), ]

  zmat      <- matrix(merged$CatNum, nrow = n_exons, ncol = n_exons, byrow = TRUE)
  hovermat  <- matrix(merged$hover,  nrow = n_exons, ncol = n_exons, byrow = TRUE)

  axis_labs <- paste0("E", ex_meta$rank)

  ## ---- 6. Plotly heatmap ----
  p <- plot_ly(
    x = axis_labs,
    y = axis_labs,
    z = zmat,
    text = hovermat,
    type = "heatmap",
    hoverinfo = "text",
    colors = cat_cols,
    zmin = 0,
    zmax = 4,
    showscale = FALSE,
    source = "mutateR_heatmap"
  )

  p <- layout(
    p,
    # REMOVED subtitle to prevent overlap
    title = if (!is.null(gene_symbol)) paste0(gene_symbol, " (", transcript_id, ")") else transcript_id,
    xaxis = list(
      title = "Exon 5'",
      categoryarray = axis_labs,
      categoryorder = "array",
      zeroline = FALSE
    ),
    yaxis = list(
      title = "Exon 3'",
      categoryarray = axis_labs,
      categoryorder = "array",
      zeroline = FALSE
    )
  )

  attr(p, "pairs_data") <- if(has_pairs) pairs_df else data.frame()

  return(p)
}
