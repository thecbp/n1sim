#' Shiny dashboard for simulating the Platform-of-1 design
#'
#' @export
runSimulator = function() {
  appDir = system.file("N1-simulator", "app", package = "n1sim")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mane`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
