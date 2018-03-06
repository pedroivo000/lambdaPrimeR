#' Open lambdaPrimeR Shiny user-interface
#'
#' @return Opens Shiny app on a new browser window
#' @export
#'
run_app_ui <- function() {
  app_dir <- system.file("lambdaPrimeRui", "shiny_app", package = 'lambdaPrimeR')
  if(app_dir == "") {
    stop("Could not find Shiny app directory. Try re-installing `lambdaPrimeR`.", call. = F)
  }
  
  shiny::runApp(app_dir, display.mode = 'normal')
}