#' GUI for Interactive Pedigree Visualization
#'
#' Running this function with no arguments will launch a general user interface for
#' visualizing pedigrees.
#'
#' @examples
#' breedr()
#'
#' @export

breedr <- function(){
  ui <- shiny::fluidPage(
    theme = bslib::bs_theme(bootswatch = "darkly"),
    shiny::titlePanel(
      shiny::h1("breedr Pedigree Dashboard", align = "center")
    ),
    shiny::fluidRow(
      shiny::column(8,
                    DiagrammeR::grVizOutput("pedigree")),
      shiny::column(4,
                    shiny::h3("Control Panel"),
                    shiny::h4("Select breedr formatted csv file"),
                    shiny::fileInput(inputId = "breedr_records",
                                     label = "Select File",
                                     placeholder = NULL),
                    shiny::conditionalPanel(condition = "output.fileUploaded",
                                            shiny::selectizeInput(inputId = "genotype",
                                                                  label = "Select Genotype",
                                                                  choices = NULL),
                                            shiny::fluidRow(
                                              shiny::column(6, shiny::selectInput(inputId = "ortn",
                                                                                  label = "Orientation",
                                                                                  choices = list('Left to Right'="LR",
                                                                                                 'Top to Bottom'="TB",
                                                                                                 'Right to Left'="RL"),
                                                                                  selected = "TB")),
                                              shiny::column(6, shiny::selectInput(inputId = "method",
                                                                                  label = "Plotting Method",
                                                                                  choices = list('Full Pedigree'="FULL",
                                                                                                 'Common Ancestry'="CA"),
                                                                                  selected = "FULL"))
                                            ),
                                            shiny::numericInput("maxn", "Max Number of Generations", 4,
                                                                min = 1, max = 100, step = 1),
                                            shiny::fluidRow(
                                              shiny::column(4, colourpicker::colourInput("col1", "Primary Color", value = "#7DFAFA")),
                                              shiny::column(4, colourpicker::colourInput("col2", "Secondary Color", value = "pink")),
                                              shiny::column(4, colourpicker::colourInput("col3", "Tertiary Color", value = "#F0B9FA"))
                                            ),
                                            shiny::downloadButton("plot_out", "Download"))
      )
    ),
    shiny::fluidRow(shiny::column(6, shiny::verbatimTextOutput("debug", T)))
  )

  server <- function(input, output, session){
    # Responds to file upload
    output$fileUploaded <- shiny::eventReactive(input$breedr_records, {
      shiny::renderText(input$breedr_records$name)
    }, ignoreInit = T)
    ds_formatted <- shiny::eventReactive(input$breedr_records, {
      read.csv(input$breedr_records$datapath)
    }, ignoreInit = T)
    geno_list <- shiny::reactive({
      ds_formatted()$Ind
    })
    shiny::observe({
      shiny::updateSelectizeInput(session, inputId = "genotype", choices = geno_list(), server = T)
    })

    # Pedigree plot reactivity
    shiny::observeEvent(input$genotype, {
      output$pedigree <- DiagrammeR::renderGrViz({
        plotigree(ds_formatted(), input$genotype, orientation = input$ortn,
                  method = input$method, n=input$maxn,
                  color1 = paste0("'",input$col1,"'"),
                  color2 = paste0("'",input$col2,"'"),
                  color3 = paste0("'",input$col3,"'"))
      })
    }, ignoreInit = T)

    # Save output with action button
    output$plot_out <- shiny::downloadHandler(
      filename = function() {
        paste0("pedigree_",
               stringr::str_replace_all(Sys.time(), " ", "_") %>%
                 stringr::str_replace_all(":", "-"),
               ".png")
      },
      content = function(file) {
        plotigree(ds_formatted(), input$genotype, orientation = input$ortn,
                  method = input$method, n=input$maxn,
                  color1 = paste0("'",input$col1,"'"),
                  color2 = paste0("'",input$col2,"'"),
                  color3 = paste0("'",input$col3,"'")) %>%
          DiagrammeRsvg::export_svg() %>% charToRaw() %>%
          rsvg::rsvg_png(file)
      })

    # Output for debugging
    output$debug <- shiny::renderPrint({
      input$col1
    })
    outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  }

  shiny::shinyApp(ui = ui, server=server)
}
