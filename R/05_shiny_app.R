
#' Shiny app for GULL
#'
#' This is the function to call a Shiny GUI app for data exploration. After preparing a GULL object, users can call this function to perform the
#' genome-wide selection and pathway specific analysis.
#'
#' @return
#' @export
#'
#' @examples
run_GULL_app <- function() {
  appDir <- system.file("shiny_app", "GULL_GUI", package = "GULL")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `GULL`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
