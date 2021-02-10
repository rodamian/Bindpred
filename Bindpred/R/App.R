#' Launches an interactive web application for labeling cells and clonotypes based on the results of an ELISA experiment, and adds the possibility to add an affintiy value for affinity prediction.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @export
#' @examples
#' \dontrun{
#' label_interactive(features = output.load_data)
#' }
#'

label_interactive <- function(features) {

    require(shiny, quietly = T)
    require(tidyverse, quietly = T)
    require(DT, quietly = T)

    features <- features %>% bind_rows %>%
        arrange(clonotype_id, decreasing = T) %>% mutate(ELISA_bind = "no") %>% distinct(cdr3s_aa, .keep_all = T)

    ui <- fluidPage(
        h1('Sequences sorted by frequency'),
        fluidRow(column(9, DT::dataTableOutput('seq_table')),
                 actionButton("save", "Save selected specific BCRs - Analyze")
        ),
        fluidRow(column(9, actionButton("close", "Close Web App and return labeled dataset"),
                           actionButton("plot", "Save selected specific BCRs - Plot repertoire features")),
                 tabPanel("Plot", fluidRow(column(5, plotOutput("AUC")),
                                           column(6, plotOutput("importance")),
                                           column(7, plotOutput("repertoire"))
                                           )
                         )
                )
    )

    server <- function(input, output) {

        output$seq_table <- renderDataTable(
            features, options = list(lengthMenu = c(5, 10, 20, 50, 100, 200),
                                     pageLength = 20,
                                     columnDefs = list(list(visible=FALSE, targets=c(10)))))

        observe({
            if(input$close > 0){
                s <- input$seq_table_rows_selected
                if (!is.null(s))
                    if (features$ELISA_bind[s] == "no")
                        features$ELISA_bind[s] <<- "yes"
                else
                    features$ELISA_bind[s] <<- "no"
                stopApp(features)
            }
        })

        observeEvent(input$plot, {
            s <- input$seq_table_rows_selected
            if (!is.null(s))
                if (features$ELISA_bind[s] == "no")
                    features$ELISA_bind[s] <<- "yes"
            else
                features$ELISA_bind[s] <<- "no"

            output$seq_table = renderDataTable(
                features, options = list(lengthMenu = c(5, 10, 20, 50, 100, 200),
                                         pageLength = 20,
                                         columnDefs = list(list(visible=FALSE, targets=c(4)))))
            output$repertoire <- renderPlot(Bindpred::explore_features(features)[[1]])
        })

        observeEvent(input$save, {
            s <- input$seq_table_rows_selected
            if (!is.null(s))
                if (features$ELISA_bind[s] == "no")
                    features$ELISA_bind[s] <<- "yes"
            else
                features$ELISA_bind[s] <<- "no"

            output$seq_table = renderDataTable(
                features, options = list(lengthMenu = c(5, 10, 20, 50, 100, 200), pageLength = 20))
            output$AUC <- renderPlot(Bindpred::classify_data(features)[[1]])
            output$importance <- renderPlot(Bindpred::classify_data(features)[[2]])
        })
    }
    runApp(shinyApp(ui = ui, server = server))
}
