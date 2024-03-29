#' @export

xplotter <- function(){
  ui <- shiny::fluidPage(
    theme = bslib::bs_theme(bootswatch = "darkly"),
    shiny::fluidRow(
      shiny::column(9,
                    shiny::h1(""),
                    ggiraph::girafeOutput(outputId = "segdist")),
      shiny::column(3,
                    shiny::h1(""),
                    shiny::fluidRow(
                      shiny::column(12, shiny::h4("Parent Options"))
                    ),
                    shiny::fluidRow(
                      shiny::column(6, shiny::selectInput(inputId = "parent1",
                                                   label = "Female",
                                                   choices = c("aa", "Aa", "AA"),
                                                   selected = "Aa")),
                      shiny::column(6, shiny::selectInput(inputId = "parent2",
                                                   label = "Male",
                                                   choices = c("aa", "Aa", "AA"),
                                                   selected = "Aa"))
                    ),
                    shiny::fluidRow(
                      shiny::column(12,
                                    shiny::h4(""),
                                    shiny::h4("Progeny Options"))
                    ),
                    shiny::fluidRow(
                      shiny::column(6, shiny::numericInput(inputId = "popsize",
                                                    label = "Population Sz.",
                                                    value = 100,
                                                    min = 1,
                                                    step = 1)),
                      shiny::column(6, shiny::numericInput(inputId = "success",
                                                    label = "Successes",
                                                    value = 1,
                                                    min = 1,
                                                    step = 1))
                    ),
                    shiny::fluidRow(
                      shiny::column(12,
                                    shiny::h4(""),
                                    shiny::h4("Model Options"))
                    ),
                    shiny::fluidRow(
                      shiny::column(6, shiny::numericInput(inputId = "ploidy",
                                                    label = "Ploidy",
                                                    value = 2,
                                                    min = 2,
                                                    max = 10,
                                                    step = 2)),
                      shiny::column(6, shiny::numericInput(inputId = "loci",
                                                    label = "Loci",
                                                    value = 1,
                                                    min = 1,
                                                    max = 4,
                                                    step = 1))
                    ), shiny::fluidRow(
                      shiny::column(12, shiny::htmlOutput("statsout"))
                    )))
  )

  server <- function(input, output, session){
    seg_df <- shiny::reactive({
      breedr::multi_seg(input$parent1, input$parent2, input$ploidy, input$loci)
    })
    ##################

    # output$debug <- shiny::renderPrint({
    #   prob()
    # })

    x <- shiny::reactive({
      ggobj = ggplot2::ggplot(seg_df() , ggplot2::aes(x=.data$Genotype, y = .data$Freq,
                                                      tooltip = paste0(Genotype, "\n", round(Freq,3)),
                                                      data_id = Genotype)) +
        ggiraph::geom_bar_interactive(stat = "identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
    })

    prob <- shiny::reactive({
      if(is.null(input$segdist_selected)){
        NULL
      } else {
        pbinom(input$success, input$popsize,
               sum(dplyr::filter(seg_df(), Genotype %in% input$segdist_selected)$Freq),
               lower.tail = F)
      }
    })

    ##################
    parent_ops <- shiny::reactive({
      multilocus_table[[input$loci]][[input$ploidy/2]]
    })
    output$statsout <- shiny::renderText({
      if(is.null(input$segdist_selected)){
        "Select genotype(s) to see binomial probabilities"
      } else {
        paste0("In a population of <b>",input$popsize,"</b> there is a <b>",
               round(prob()*100, 3),"%</b> chance of recovering at least <b>",input$success,
               "</b> of the following genotype(s): <b>", paste0(input$segdist_selected, collapse = ", "),"</b>")
      }
    })
    shiny::observeEvent(list(input$ploidy, input$loci), {
      shiny::updateSelectInput(session, "parent1",
                               choices = parent_ops(),
                               selected = parent_ops()[ceiling(length(parent_ops())/2)])
      shiny::updateSelectInput(session, "parent2",
                               choices = parent_ops(),
                               selected = parent_ops()[ceiling(length(parent_ops())/2)])
    }, ignoreInit = T)
    output$segdist <- ggiraph::renderGirafe({
      ggiraph::girafe(ggobj = x(),
                      options = list(ggiraph::opts_hover(css = "fill:#00FFFF;stroke:gray;r:5pt;",
                                                         reactive = T),
                                     ggiraph::opts_selection(css = "fill:#0096FF;")),
                      width_svg = 9
      )
    })
  }
  shiny::shinyApp(ui = ui, server=server)
}
