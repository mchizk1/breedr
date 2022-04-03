
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
                                            shiny::column(6, shinyWidgets::radioGroupButtons(inputId = "ortn",
                                                                                             label = "Orientation",
                                                                                             choices = c("LR", "TB", "RL"),
                                                                                             selected = "TB")),
                                            shiny::column(6, shinyWidgets::radioGroupButtons(inputId = "method",
                                                                                             label = "Plotting Method",
                                                                                             choices = c("FULL", "CA"),
                                                                                             selected = "FULL"))
                                          ),
                                          shiny::numericInput("maxn", "Max number of Generations", 4,
                                                              min = 1, max = 100, step = 1),
                                          shiny::actionButton("saveplot", "Save Pedigree"))
                  )
  ),
  shiny::fluidRow(shiny::column(6, shiny::verbatimTextOutput("debug", T)))
)
