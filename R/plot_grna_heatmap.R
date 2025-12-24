#' @title Plot exon‑phase compatibility as a discrete heatmap
#'
#' @description Helper function for large-gene visualisation.
#'
#' @param exon_gr GRanges from get_exon_structures(output="GRanges")
#' @param transcript_id Character. Transcript label.
#' @param gene_symbol Character. Gene label.
#' @param pairs_df Optional data.frame of gRNA‑pair information
#'        (currently unused but accepted for forward compatibility).
#'
#' @return ggplot object
#' @export
plot_grna_heatmap <- function(exon_gr,
                              transcript_id = "Transcript",
                              gene_symbol = NULL,
                              species,
                              pairs_df = NULL) {

  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(RColorBrewer)
  })

  ## ----- 1. Extract exon metadata -----
  ex_meta <- as.data.frame(mcols(exon_gr))
  ex_meta$rank <- seq_len(nrow(ex_meta))
  ex_meta <- ex_meta %>%
    mutate(is_UTR = (start_phase == -1 & end_phase == -1))

  ## ----- 2. Compute phase compatibility -----
  comp_df <- check_exon_phase(ex_meta, include_contiguous = TRUE)
  comp_df$Category <- ifelse(comp_df$compatible, "Compatible", "Incompatible")

  all_pairs <- expand.grid(exon_5p = ex_meta$rank,
                           exon_3p = ex_meta$rank)
  mat_df <- left_join(all_pairs, comp_df,
                      by = c("exon_5p", "exon_3p")) %>%
    mutate(Category = case_when(
      exon_5p == exon_3p                ~ "Self",
      exon_3p - exon_5p == 1            ~ "Contiguous",
      is.na(Category)                   ~ "Incompatible",
      TRUE                              ~ Category
    ))

  utr_exons <- ex_meta %>% filter(is_UTR) %>% pull(rank)
  mat_df <- mat_df %>%
    mutate(Category = ifelse(exon_5p %in% utr_exons |
                               exon_3p %in% utr_exons,
                             "UTR", Category))

  mat_df$Category <- factor(mat_df$Category,
                            levels = c("Self","Contiguous",
                                       "Compatible","Incompatible","UTR"))

  cat_cols <- c(
    "Self" = "grey70",
    "Contiguous" = "palegreen3",
    "Compatible"  = "goldenrod2",
    "Incompatible" = "steelblue4",
    "UTR" = "firebrick3"
  )

  ## ----- 3. Retrieve domain annotations -----
  domain_df <- tryCatch(
    map_protein_domains(transcript_id,
                        species = species,
                        source = "pfam"),
    error = function(e) {
      message("Domain retrieval failed for transcript ", transcript_id, ": ", e$message)
      NULL
    }
  )

  domain_plot_df <- data.frame()
  if (!is.null(domain_df) && nrow(domain_df) > 0) {
    domain_plot_df <- domain_df %>%
      group_by(domain_desc) %>%
      summarise(exon_start = min(exon_rank, na.rm = TRUE),
                exon_end   = max(exon_rank, na.rm = TRUE),
                .groups = "drop") %>%
      distinct(domain_desc, .keep_all = TRUE)
  }

  ## ----- 4. Base heatmap -----
  p_heat <- ggplot(mat_df, aes(x = exon_5p, y = exon_3p, fill = Category)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = cat_cols,
                      name = "Phase\ncompatibility") +
    coord_fixed() +
    scale_x_continuous(breaks = ex_meta$rank,
                       labels = paste0("E", ex_meta$rank),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = ex_meta$rank,
                       labels = paste0("E", ex_meta$rank),
                       expand = c(0, 0)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.margin = margin(6, 6, 25, 6),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(paste0(
      if (!is.null(gene_symbol))
        paste0(transcript_id, " (", gene_symbol, ")")
      else transcript_id,
      " - Exon phase-compatibility matrix"
    ))

  ## ----- 5. Domain annotation layer -----
  if (!is.null(domain_plot_df) && nrow(domain_plot_df) > 0) {
    n_domains <- nrow(domain_plot_df)

    # Vertical tracks (all below the heatmap)
    domain_plot_df <- domain_plot_df %>%
      mutate(
        y_base = -1.2,
        offset = seq_len(n_domains) * 0.35,
        ymin = y_base - offset,
        ymax = y_base - offset + 0.25
      )

    domain_cols <- RColorBrewer::brewer.pal(min(max(3, n_domains), 8), "Set2")

    # Assign label placement: alternate above/below
    domain_plot_df$label_y <- ifelse(seq_len(n_domains) %% 2 == 1,
                                     domain_plot_df$ymax + 0.1,  # odd -> above
                                     domain_plot_df$ymin - 0.1)  # even -> below

    # ---- Added safety guard ----
    if (nrow(domain_plot_df) > 0 && length(domain_cols) > 0) {
      p_heat <- p_heat +
        geom_rect(
          data = domain_plot_df,
          aes(xmin = exon_start - 0.45,
              xmax = exon_end + 0.45,
              ymin = ymin,
              ymax = ymax),
          inherit.aes = FALSE,
          fill = domain_cols[seq_len(min(length(domain_cols), n_domains))],
          colour = "black",
          linewidth = 0.25,
          alpha = 0.85
        ) +
        geom_text(
          data = domain_plot_df,
          aes(x = (exon_start + exon_end)/2,
              y = label_y,
              label = domain_desc),
          inherit.aes = FALSE,
          size = 3.1,
          fontface = "italic",
          colour = "black"
        ) +
        expand_limits(y = min(domain_plot_df$ymin) - 0.3)
    }
  } else {
    message("No protein domain annotations available; displaying heatmap only.")
  }

  return(p_heat)
}
