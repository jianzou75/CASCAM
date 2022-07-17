
#' Shiny app for CASCAM
#'
#' This is the function to call a Shiny GUI app for data exploration. After preparing a CASCAM object, users can call this function to perform the
#' genome-wide selection and pathway specific analysis.
#'
#' @return
#' @export
#'
#' @examples
run_CASCAM_app <- function() {
  appDir <- system.file("shiny_app", "CASCAM_GUI", package = "CASCAM")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `CASCAM`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
