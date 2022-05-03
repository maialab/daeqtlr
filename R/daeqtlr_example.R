#' Get path to daeqtlr example
#'
#' daeqtlr comes bundled with a number of sample files in its `inst/extdata`
#' directory. This function make them easy to access.
#'
#' @param file Name of file. If `NULL`, the example files will be listed.
#' @export
#'
#' @md
#' @examples
#' daeqtlr_example()
#' daeqtlr_example("zygosity.csv")
daeqtlr_example <- function(file = NULL) {
  if (is.null(file)) {
    dir(system.file("extdata", package = "daeqtlr"))
  } else {
    system.file("extdata", file, package = "daeqtlr", mustWork = TRUE)
  }
}
