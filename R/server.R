
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
                method = input$method, n=input$maxn)
    })
  }, ignoreInit = T)

  # Save output with action button
  shiny::observeEvent(input$saveplot, {
    file.choose()
  })

  # Output for debugging
  output$debug <- shiny::renderPrint({
    input$ortn
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
}
