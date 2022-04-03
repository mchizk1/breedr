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
  shiny::shinyApp(ui = ui, server=server)
}
