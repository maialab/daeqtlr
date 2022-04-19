#' Import DAE SNPs
#'
#' This function reads in a list of DAE SNPs from a text file. One SNP per row.
#' By default, it expects a header, which will be ignored anyway. If the first
#' line is immediately a SNP, then use \code{header = FALSE}.
#'
#' @param file A path to a file.
#' @param header Is the first line a line with a column name?
#' @param distinct Remove any duplicate SNPs?
#'
#' @return A one-column data frame. This single column is named \code{dae_snp}.
#'   Rows are the DAE SNPs read from the text file.
#'
#' @export
read_dae_snps <- function(file, header = TRUE, distinct = TRUE) {
  read_snps(file, colname = 'dae_snp', header = header, distinct = distinct)
}
