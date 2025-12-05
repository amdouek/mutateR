#' Interactive viewer for mutateR gRNA design heatmaps
#'
#' Launches a small Shiny gadget showing the interactive Plotly heatmap.
#' Clicking a cell filters the gRNA table below. Primer3-generated genotyping primers are displayed for each recommended gRNA pair.
#'
#' @param plot_obj Plotly object returned in the mutateR output list
#'                 (element \code{$plot}); must have attribute \code{"pairs_data"}.
#' @export
mutateR_viewer <- function(plot_obj) {

  if (!requireNamespace("shiny", quietly = TRUE) ||
      !requireNamespace("plotly", quietly = TRUE) ||
      !requireNamespace("DT", quietly = TRUE)) {
    stop("Packages 'shiny', 'plotly', and 'DT' are required for mutateR_viewer().")
  }

  # ----- Guard clause: Ensure input is interactive -----
  if (inherits(plot_obj, "ggplot")) {
    stop("mutateR_viewer() requires an interactive plotly object.", call. = FALSE)
  }

  if (!inherits(plot_obj, "plotly")) {
    stop("Input must be a plotly object.", call. = FALSE)
  }

  pairs_df <- attr(plot_obj, "pairs_data")
  # Handle empty case gracefully
  if (is.null(pairs_df)) pairs_df <- data.frame()

  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::tags$head(shiny::tags$style(
        shiny::HTML(".shiny-output-error-validation { color: red; }")
      )),
      shiny::titlePanel("mutateR - Interactive Design Viewer"),

      shiny::sidebarLayout(
        shiny::sidebarPanel(
          width = 3,
          shiny::h4("Selection"),
          shiny::htmlOutput("sel_summary"),
          shiny::hr(),
          shiny::helpText("1. Click a cell in the heatmap."),
          shiny::helpText("2. Review candidate pairs."),
          shiny::helpText("3. Select rows in table to export."),
          shiny::hr(),

          # --- Export Section ---
          shiny::h4("Data Export"),
          shiny::downloadButton("dl_csv", "Download Table (CSV)"),
          shiny::br(), shiny::br(),

          shiny::hr(),
          shiny::h4("Vendor Export"),
          shiny::selectInput("vendor_format", "Format:",
                             choices = c("IDT (Bulk Upload)", "Generic (Name, Seq)")),

          shiny::helpText("Exports gRNA protospacer sequences and genotyping primers."),
          shiny::downloadButton("dl_vendor", "Download Order Sheet"),

          shiny::hr(),
          shiny::actionButton("close_btn", "Close Viewer", class = "btn-danger")
        ),
        shiny::mainPanel(
          width = 9,
          plotly::plotlyOutput("heatmap", height = "600px"),
          shiny::hr(),
          shiny::h4("Candidate gRNA Pairs"),
          DT::dataTableOutput("pairs_tbl")
        )
      )
    ),

    server = function(input, output, session) {

      # Store filtered data based on heatmap click
      heatmap_filtered_data <- shiny::reactiveVal(data.frame())

      output$heatmap <- plotly::renderPlotly({
        plotly::event_register(plot_obj, "plotly_click")
        plot_obj
      })

      # --- Heatmap Click Event ---
      shiny::observeEvent(plotly::event_data("plotly_click", source = "mutateR_heatmap"), {
        click <- plotly::event_data("plotly_click", source = "mutateR_heatmap")
        if (is.null(click)) return()

        x_lab <- click$x
        y_lab <- click$y
        e5 <- suppressWarnings(as.integer(gsub("E", "", x_lab)))
        e3 <- suppressWarnings(as.integer(gsub("E", "", y_lab)))

        if (is.na(e5) || is.na(e3)) return()

        if (nrow(pairs_df) > 0) {
          sub_df <- pairs_df %>%
            dplyr::filter(upstream_pair == e5, downstream_pair == e3) %>%
            # Sort by recommended, then average score
            dplyr::arrange(dplyr::desc(recommended),
                           dplyr::desc((ontarget_score_5p + ontarget_score_3p)/2))
          heatmap_filtered_data(sub_df)
        } else {
          heatmap_filtered_data(data.frame())
        }

        output$sel_summary <- shiny::renderUI({
          dat <- heatmap_filtered_data()
          if(nrow(dat) > 0) {
            n_total <- nrow(dat)
            n_rec   <- sum(dat$recommended, na.rm = TRUE)
            shiny::HTML(paste0(
              "<b>Exon Pair:</b> E", e5, " - E", e3, "<br/>",
              "<b>Candidate Pairs:</b> ", n_total, " (", n_rec, " recommended)"
            ))
          } else {
            shiny::HTML(paste0(
              "<b>Exon Pair:</b> E", e5, " - E", e3, "<br/>",
              "<b>Candidate Pairs:</b> 0"
            ))
          }
        })
      })

      # --- Data Table ---
      output$pairs_tbl <- DT::renderDataTable({
        dat <- heatmap_filtered_data()

        if (nrow(dat) == 0) {
          return(DT::datatable(data.frame(Info = "Click a heatmap cell containing data to view gRNAs."),
                               options = list(dom = 't'), rownames = FALSE))
        }

        # Select relevant columns for display
        cols_to_show <- c("exon_5p", "exon_3p", "genomic_deletion_size",
                          "protospacer_sequence_5p", "ontarget_score_5p",
                          "protospacer_sequence_3p", "ontarget_score_3p",
                          "recommended") # For consideration: use the 'colvis' button option to allow users to select/deselect other columns from the $pairs dataframe

        # Add primer info if it exists
        if("primer_ext_fwd" %in% names(dat)) {
          cols_to_show <- c(cols_to_show, "primer_ext_fwd", "primer_ext_rev", "primer_int_fwd", "primer_int_rev")
        }

        DT::datatable(dat[, intersect(names(dat), cols_to_show), drop=FALSE],
                      extensions = c("Buttons", "Select"),
                      rownames = FALSE,
                      selection = 'multiple', # Allow multiple row selection. Works but shiny doesn't like DT's selection implementation - for consideration
                      options = list(
                        dom = 'Blfrtip',
                        buttons = list('copy','csv'),
                        lengthMenu = list(c(5, 10, 25), c('5', '10', '25')),
                        scrollX = TRUE,
                        select = list(style = 'multi')
                      )) %>%
          DT::formatRound(columns = intersect(names(dat), c("ontarget_score_5p", "ontarget_score_3p")),
                          digits = 2)
      }, server = FALSE)

      # --- CSV data download ---
      output$dl_csv <- shiny::downloadHandler(
        filename = function() { paste0("mutateR_data_", Sys.Date(), ".csv") },
        content = function(file) {
          # Download rows currently in the table (filtered by heatmap click)
          utils::write.csv(heatmap_filtered_data(), file, row.names = FALSE)
        }
      )

      # --- Vendor export logic ---
      output$dl_vendor <- shiny::downloadHandler(
        filename = function() {
          paste0("mutateR_IDT_Order_", Sys.Date(), ".csv")
        },
        content = function(file) {
          # 1. Get data: Use rows SELECTED in the table. If none selected, use ALL in table. <- Consider later if this is the best course of action... after the backend filtering we've applied this shouldn't ever be more than an easily correctable nuisance for users who error (minimal crash/memory fuckery risk unlikely for a subset of a subset) BUT maybe explicitly declare default behaviour or switch to 'no rows selected' error and cognate behaviour.
          full_data <- heatmap_filtered_data()
          s <- input$pairs_tbl_rows_selected

          if (length(s) > 0) {
            export_data <- full_data[s, , drop=FALSE]
          } else {
            export_data <- full_data
          }

          if (nrow(export_data) == 0) return(NULL)

          # 2. Build order list (currently separates oligos and guides by flagging with a DNA/RNA tag, but this is not suitable for the vendor bulk upload and would need to go)
          order_rows <- list()

          for(i in seq_len(nrow(export_data))) {
            row <- export_data[i, ]
            # Append sequence name (mainly for primers - need to figure out an elegant way to handle crRNA seq export as this is distinct from ordering oligos). Also need to factor in common name length restrictions (i.e. make a standardised oligo nomenclature -- think on this)
            id_base <- paste0("Pair_", i, "_E", row$exon_5p, "-E", row$exon_3p)

            # --- 5' gRNA (RNA) ---
            order_rows[[length(order_rows)+1]] <- data.frame(
              Name = paste0(id_base, "_gRNA_5p"),
              Sequence = row$protospacer_sequence_5p,
              Molecule = "RNA (crRNA)",
              Scale = "2nm", Purification = "STD" # Note that we are automatically applying scale and purification method, which might be risky (don't want the user to accidentally order the wrong amount etc.). Should make this a user choice, but I don't want to make the interface too busy with dropdown boxes etc. etc.)
            )

            # --- 3' gRNA (RNA) ---
            order_rows[[length(order_rows)+1]] <- data.frame(
              Name = paste0(id_base, "_gRNA_3p"),
              Sequence = row$protospacer_sequence_3p,
              Molecule = "RNA (crRNA)",
              Scale = "2nm", Purification = "STD"
            )

            # --- External Primers (DNA) -- needed for all gRNA pairs ---
            if ("primer_ext_fwd" %in% names(row) && !is.na(row$primer_ext_fwd)) {
              order_rows[[length(order_rows)+1]] <- data.frame(
                Name = paste0(id_base, "_Prim_Ext_F"), Sequence = row$primer_ext_fwd,
                Molecule = "DNA (Oligo)", Scale = "25nm", Purification = "STD"
              )
              order_rows[[length(order_rows)+1]] <- data.frame(
                Name = paste0(id_base, "_Prim_Ext_R"), Sequence = row$primer_ext_rev,
                Molecule = "DNA (Oligo)", Scale = "25nm", Purification = "STD"
              )
            }

            # --- Internal Primers (DNA) - Strategy B (only for very large deletions) ---
            if ("primer_int_fwd" %in% names(row) && !is.na(row$primer_int_fwd)) {
              order_rows[[length(order_rows)+1]] <- data.frame(
                Name = paste0(id_base, "_Prim_Int_F"), Sequence = row$primer_int_fwd,
                Molecule = "DNA (Oligo)", Scale = "25nm", Purification = "STD"
              )
              order_rows[[length(order_rows)+1]] <- data.frame(
                Name = paste0(id_base, "_Prim_Int_R"), Sequence = row$primer_int_rev,
                Molecule = "DNA (Oligo)", Scale = "25nm", Purification = "STD"
              )
            }
          }

          final_df <- do.call(rbind, order_rows)

          # 3. Write output
          utils::write.csv(final_df, file, row.names = FALSE)
        }
      )

      shiny::observeEvent(input$close_btn, {
        shiny::stopApp()
      })
    }
  )
}
