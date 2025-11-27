#' Plot exon structure, phase compatibility and recommended gRNA pairs
#'
#' Displays three schematic layers:
#'   • Top dashed grey arcs – phase‑compatible exon pairs
#'   • Middle blue rectangles – exons
#'   • Bottom coloured arcs – recommended gRNA pair targets
#' Left‑side labels identify each layer.
#'
#' @param exon_gr       GRanges from get_exon_structures(output="GRanges")
#' @param pairs_df      Data frame from assemble_grna_pairs()
#' @param transcript_id Character transcript/gene label for the title.
#' @param top_n         Integer number of recommended pairs to display - arc mode only.
#'                      Default NULL (shows all recommended pairs).
#' @param mode          Character, one of "heat" (default) or "arc".
#' @param interactive Logical; default FALSE. If TRUE, launches an interactive plotly-based viewer allowing hover-over display of exon-pair and gRNA-pair data, multiple selection, and CSV export.
#'
#' @return A ggplot object
#' @export
plot_grna_design <- function(exon_gr,
                             pairs_df = NULL,
                             transcript_id = "Transcript",
                             gene_symbol = NULL,
                             species,
                             top_n = NULL,
                             mode = c('heat', 'arc'),
                             interactive = FALSE) {
  library(ggplot2); library(dplyr); library(tibble); library(RColorBrewer)

  mode <- match.arg(mode)

  ## ---- Define heatmap mode -----
  if (mode == 'heat') {
    message("Plotting exon phase compatibility and gRNA pairs...")
    if (!exists("plot_grna_heatmap", mode = "function")) {
      stop("The function plot_grna_heatmap() must be defined to use heatmap mode.\n Ensure the package has been properly installed.")
    }
    # Pass pairs_df for gRNA-density overlay -- placeholder for future development
    return(
      plot_grna_heatmap(exon_gr = exon_gr,
                        transcript_id = transcript_id,
                        gene_symbol = gene_symbol,
                        species = species,
                        pairs_df = pairs_df)
    )
  }

  ## ---- Define arc mode -----
  message("Plotting exon phase compatibility and gRNA pairs...")

  exon_df <- as.data.frame(exon_gr)
  n_exons <- nrow(exon_df)
  exon_df$rank <- seq_len(n_exons)

  ## ---- Early exit for large genes --------------------------------------
  if (length(exon_gr) > 12) {
    message("Detected large gene (>12 exons); switching to heatmap mode.")
    return(
      plot_grna_heatmap(exon_gr = exon_gr,
                        transcript_id = transcript_id,
                        gene_symbol = gene_symbol,
                        pairs_df = pairs_df)
      )
  }

  ## -- Detect intragenic mode ------------------------------------------
  intragenic_mode <- n_exons <= 2 ||
    (!is.null(pairs_df) &&
       "del_size" %in% names(pairs_df) &&
       all(pairs_df$exon_5p == pairs_df$exon_3p |
             abs(pairs_df$exon_3p - pairs_df$exon_5p) <= 1))

  ## ----- 1. Intragenic plotting mode (for genes with (≤2 exons)) -----
  if (intragenic_mode) {
    message("Plotting intragenic deletion mode.")
    library(ggplot2); library(dplyr); library(tibble)
    # Build base exon box(es)
    exon_blocks <- tibble(
      x = seq_len(n_exons),
      xmin = seq_len(n_exons) - 0.4,
      xmax = seq_len(n_exons) + 0.4,
      ymin = 0, ymax = 0.5,
      label = paste0("E", seq_len(n_exons))
    )

    # Select pairs to show
    if (is.null(pairs_df) || nrow(pairs_df) == 0)
      pairs_df <- data.frame()
    if (!is.null(top_n) && nrow(pairs_df) > top_n)
      pairs_df <- head(pairs_df, top_n)

    # Define colour palette
    n_pairs <- nrow(pairs_df)
    cols <- brewer.pal(min(max(3, n_pairs), 12), "Set3")
    if (n_pairs > 12) cols <- rep(cols, length.out = n_pairs)

    # Define each deletion span as an arc or diagonal line
    if (n_pairs > 0) {
      arc_tbl <- bind_rows(lapply(seq_len(n_pairs), function(i) {
        df <- pairs_df[i, ]
        max_del <- max(pairs_df$del_size, na.rm = TRUE)
        rel <- if (!is.na(df$del_size) && max_del > 0)
          0.3 + 0.7 * (df$del_size / max_del) else 0.3
        tibble(
          pair_id = paste0("Pair_", i),
          x = 1 + 0.1,
          xend = n_exons - 0.1,
          y = 0.5 + rel * 0.3,
          yend = 0.5 + rel * 0.3,
          label = paste0("Δ", round(df$del_size / 1000, 1), " kb")
        )
      }))
      arc_tbl$col <- cols[seq_len(nrow(arc_tbl))]
    } else {
      arc_tbl <- tibble()
    }

    # --- then plotting ---
    p <- ggplot() +
      geom_rect(data = exon_blocks,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "steelblue3", colour = "black") +
      geom_text(data = exon_blocks,
                aes(x = x, y = ymax + 0.05, label = label),
                size = 5)

    if (nrow(arc_tbl) > 0) {
      p <- p +
        geom_segment(data = arc_tbl,
                     aes(x = x, xend = xend, y = y, yend = yend, colour = pair_id),
                     linewidth = 1.5, lineend = "round") +
        geom_text(data = arc_tbl,
                  aes(x = (x + xend)/2, y = y + 0.05, label = label),
                  size = 4, colour = "black")
    }

    title_txt <- paste0(
      if (!is.null(gene_symbol))
        paste0(gene_symbol, " (", transcript_id, ")")
      else transcript_id,
      " — Intragenic large-deletion design (≤2 exons)"
    )

    p <- p +
      scale_colour_manual(values = cols) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      coord_cartesian(clip = "off") +
      ggtitle(title_txt)

    return(p)
  }


  ## ----- 2. Multi-exon arc mode -----
  library(ggplot2); library(dplyr); library(tibble); library(RColorBrewer)

  no_pairs <- is.null(pairs_df) || nrow(pairs_df) == 0
  if (no_pairs) {
    message("No recommended gRNA pairs available. Plotting exon structure and phase-compatibility only.")
    pairs_df <- data.frame()
  }

  exon_df <- exon_df %>%
    arrange(rank) %>%
    mutate(label = paste0("E", rank),
           xmin = rank - 0.45, xmax = rank + 0.45,
           ymin = 0, ymax = 0.5)

  make_arc <- function(x1, x2, h = 0.3, n = 50, offset = 0) {
    tibble(
      x = seq(x1, x2, length.out = n),
      y = h * sin(seq(0, pi, length.out = n)) + offset
    )
  }

  comp_pairs <- data.frame()
  if (length(exon_gr) > 1) {
    comp_pairs <- check_exon_phase(as.data.frame(mcols(exon_gr)),
                                   include_contiguous = FALSE)
    comp_pairs <- subset(comp_pairs, compatible == TRUE)
  }
  if (nrow(comp_pairs) > 0) {
    top_arcs <- bind_rows(
      lapply(seq_len(nrow(comp_pairs)), function(i) {
        e5 <- comp_pairs$exon_5p[i]; e3 <- comp_pairs$exon_3p[i]
        make_arc(e5, e3, h = 0.25 + 0.05 * i, offset = 0.6) %>%
          mutate(label = paste0("E", e5, "–E", e3))
      })
    )
  } else {
    top_arcs <- tibble(x = numeric(), y = numeric(), label = character())
  }

  bottom_arcs <- tibble(x = numeric(), y = numeric(), label = character(), pair_id = character())
  pair_cols <- c()
  if (!no_pairs) {
    pairs_df <- pairs_df %>%
      filter(recommended == TRUE,
             !is.na(exon_5p), !is.na(exon_3p),
             !is.na(upstream_pair), !is.na(downstream_pair)) %>%
      arrange(upstream_pair, downstream_pair) %>%
      mutate(pair_id = paste0("Pair_", seq_len(n())))
    if (!is.null(top_n)) pairs_df <- slice_head(pairs_df, n = top_n)

    n_pairs <- nrow(pairs_df)
    base_cols <- brewer.pal(min(max(3, n_pairs), 12), "Set3")
    if (n_pairs > 12) base_cols <- rep(base_cols, length.out = n_pairs)
    pair_cols <- setNames(base_cols, pairs_df$pair_id)

    bottom_arcs <- bind_rows(lapply(seq_len(n_pairs), function(i) {
      e5 <- pairs_df$exon_5p[i]; e3 <- pairs_df$exon_3p[i]
      make_arc(e5, e3, h = -0.25, offset = -0.25) %>%
        mutate(pair_id = pairs_df$pair_id[i],
               label = paste0("E", e5, "+E", e3))
    }))
  }

  label_df <- tibble(
    x = 0.2,
    y = c(max(top_arcs$y, na.rm = TRUE) - 0.5,
          0.2,
          if (no_pairs) -0.2 else min(bottom_arcs$y, na.rm = TRUE) - 0.05),
    label = c("Phase-compatible exons", "Exons",
              if (no_pairs) "No recommended pairs" else "Recommended targets")
  )

  p <- ggplot(exon_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "steelblue3", colour = "black") +
    geom_text(aes(x = rank, y = ymax + 0.05, label = label), size = 5) +
    geom_path(data = top_arcs, aes(x = x, y = y, group = label),
              colour = "grey50", linetype = "dashed", linewidth = 1) +
    geom_text(data = top_arcs %>% group_by(label) %>%
                summarise(x = mean(range(x)), y = max(y)),
              aes(x = x, y = y + 0.05, label = label),
              size = 5, colour = "grey30")

  if (!no_pairs && nrow(bottom_arcs) > 0) {
    p <- p +
      geom_path(data = bottom_arcs,
                aes(x = x, y = y, group = label, colour = pair_id),
                linewidth = 1.4) +
      geom_text(data = bottom_arcs %>%
                  group_by(label) %>%
                  summarise(x = mean(range(x)), y = min(y)),
                aes(x = x, y = y - 0.05, label = label),
                size = 5, colour = "black") +
      scale_colour_manual(values = pair_cols)
  }

  p <- p +
    geom_text(data = label_df,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.5,
              size = 4.2, fontface = "bold", colour = "black") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(paste0(
      if (!is.null(gene_symbol))
        paste0(gene_symbol, " (", transcript_id, ")")
      else transcript_id,
      " — Phase-compatibility", if (no_pairs) "" else " & recommended gRNA pairs"
    )) +
    coord_cartesian(clip = "off")

  return(p)
}
