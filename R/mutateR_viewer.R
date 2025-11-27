#' Interactive viewer for mutateR gRNA design heatmaps
#'
#' Launches a small Shiny gadget showing the interactive Plotly heatmap.
#' Clicking a cell filters the gRNA table below.
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
    stop("mutateR_viewer() requires an interactive plotly object.\n",
         "The provided object is a static ggplot (the default output of run_mutateR).\n",
         "Please generate an interactive plot (e.g., using plot_grna_interactive(), or specifying interactive = TRUE in run_mutateR) before passing it to this viewer.",
         call. = FALSE)
  }

  if (!inherits(plot_obj, "plotly")) {
    stop("Input must be a 'plotly' object.", call. = FALSE)
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
          # Changed to htmlOutput to support line breaks
          shiny::htmlOutput("sel_summary"),
          shiny::hr(),
          shiny::helpText("1. Click a cell in the heatmap."),
          shiny::helpText("2. Review candidate pairs below."),
          shiny::helpText("3. Select rows to export specific pairs."),
          shiny::hr(),
          shiny::h5("Export gRNAs from selected exon pair"),
          shiny::downloadButton("dl_csv", "Download all gRNA pairs"),
          shiny::br(), shiny::br(),
          shiny::downloadButton("dl_rec", "Download recommended gRNA pairs"),
          shiny::br(), shiny::br(),
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

      selected_data <- shiny::reactiveVal(data.frame())

      output$heatmap <- plotly::renderPlotly({
        plotly::event_register(plot_obj, "plotly_click")
        plot_obj
      })

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
            dplyr::arrange(dplyr::desc((ontarget_score_5p + ontarget_score_3p)/2))
          selected_data(sub_df)
        } else {
          selected_data(data.frame())
        }

        # Changed to renderUI to support HTML tags
        output$sel_summary <- shiny::renderUI({
          dat <- selected_data()
          if(nrow(dat) > 0) {
            n_total <- nrow(dat)
            n_rec   <- sum(dat$recommended, na.rm = TRUE)
            shiny::HTML(paste0(
              "<b>Exon Pair:</b> E", e5, " – E", e3, "<br/>",
              "<b>Candidate gRNA Pairs:</b> ", n_total, " (", n_rec, " recommended)"
            ))
          } else {
            shiny::HTML(paste0(
              "<b>Exon Pair:</b> E", e5, " – E", e3, "<br/>",
              "<b>Candidate gRNA Pairs:</b> 0"
            ))
          }
        })
      })

      output$pairs_tbl <- DT::renderDataTable({
        dat <- selected_data()

        if (nrow(dat) == 0) {
          return(DT::datatable(data.frame(Info = "Click a heatmap cell containing data to view gRNAs."),
                               options = list(dom = 't'), rownames = FALSE))
        }

        desired_cols <- c("exon_5p", "exon_3p",
                          "genomic_deletion_size", "transcript_deletion_size",
                          "protospacer_sequence_5p", "ontarget_score_5p",
                          "protospacer_sequence_3p", "ontarget_score_3p",
                          "recommended")

        show_cols <- intersect(names(dat), desired_cols)

        DT::datatable(dat[, show_cols, drop=FALSE],
                      extensions = c("Buttons", "Select"),
                      rownames = FALSE,
                      selection = 'none',
                      options = list(
                        dom = 'Blfrtip',
                        lengthMenu = list(c(10, 25, 50, 100, -1),
                                          c('10', '25', '50', '100', 'All')),
                        scrollX = TRUE,
                        select = list(style = 'multi'),
                        buttons = list(
                          list(extend = 'copy',
                               text = 'Copy Selected',
                               exportOptions = list(modifier = list(selected = TRUE))),
                          list(extend = 'csv',
                               text = 'Export Selected (CSV)',
                               exportOptions = list(modifier = list(selected = TRUE)))
                        )
                      )) %>%
          DT::formatRound(columns = intersect(names(dat), c("ontarget_score_5p", "ontarget_score_3p")),
                          digits = 3)
      }, server = FALSE)

      output$dl_csv <- shiny::downloadHandler(
        filename = function() { paste0("mutateR_all_", Sys.Date(), ".csv") },
        content = function(file) {
          utils::write.csv(selected_data(), file, row.names = FALSE)
        }
      )

      output$dl_rec <- shiny::downloadHandler(
        filename = function() { paste0("mutateR_recommended_", Sys.Date(), ".csv") },
        content = function(file) {
          dat <- selected_data()
          if (nrow(dat) > 0 && "recommended" %in% names(dat)) {
            dat <- dat[dat$recommended == TRUE, ]
          }
          utils::write.csv(dat, file, row.names = FALSE)
        }
      )

      shiny::observeEvent(input$close_btn, {
        shiny::stopApp()
      })
    }
  )
}
