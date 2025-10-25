#' Plot exon‑phase compatibility as a discrete heatmap
#' (helper function for large-gene visualisation)
#'
#' @param exon_gr GRanges from get_exon_structures(output="GRanges")
#' @param transcript_id Character. Transcript label.
#' @param gene_symbol Character. Gene label.
#'
#' @return ggplot object
#' @export
plot_grna_heatmap <- function(exon_gr, 
                              transcript_id = "Transcript",
                              gene_symbol = NULL) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })
  
  # ----- 1. Extract exon metadata and flag UTR exons -----
  ex_meta <- as.data.frame(mcols(exon_gr))
  ex_meta$rank <- seq_len(nrow(ex_meta))
  ex_meta <- ex_meta %>%
    mutate(is_UTR = (start_phase == -1 & end_phase == -1))
  
  # ----- 2. Compute phase compatibility between exons -----
  comp_df <- check_exon_phase(ex_meta, include_contiguous = TRUE)
  comp_df$Category <- ifelse(comp_df$compatible, "Compatible", "Incompatible")
  
  # ----- 3. Build full exon pair matrix and classify cells -----
  all_pairs <- expand.grid(exon_5p = ex_meta$rank,
                           exon_3p = ex_meta$rank)
  
  mat_df <- left_join(all_pairs, comp_df,
                      by = c("exon_5p", "exon_3p")) %>%
    mutate(Category = case_when(
      exon_5p == exon_3p                ~ "Self",
      exon_3p - exon_5p == 1            ~ "Contiguous",  # new category
      is.na(Category)                   ~ "Incompatible",
      TRUE                              ~ Category
    ))
  
  # Flag any pairs involving UTR exons
  utr_exons <- ex_meta %>% filter(is_UTR) %>% pull(rank)
  mat_df <- mat_df %>%
    mutate(Category = ifelse(exon_5p %in% utr_exons |
                               exon_3p %in% utr_exons,
                             "UTR", Category))
  
  # ------ 4. Ensure factor levels & colour mapping -----
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
  
  # ----- 5. Assemble heatmap -----
  ggplot(mat_df, aes(x = exon_5p, y = exon_3p, fill = Category)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = cat_cols,
                      name = "Phase\ncompatibility") +
    coord_fixed() +
    scale_x_continuous(breaks = ex_meta$rank,
                       labels = paste0("E", ex_meta$rank)) +
    scale_y_continuous(breaks = ex_meta$rank,
                       labels = paste0("E", ex_meta$rank)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(paste0(
      if (!is.null(gene_symbol))
        paste0(transcript_id, " (", gene_symbol, ")")
      else transcript_id,
      " — Exon phase-compatibility matrix"
    ))
}