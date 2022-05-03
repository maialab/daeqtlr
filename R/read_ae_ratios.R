#' Import expression ratios
#'
#' This function reads in expression ratios from a text file. The data is
#' expected to be in tabular format. First column is expected to be the locus
#' identity (e.g. SNP identifier), and remaining columns are the samples.
#'
#' @param file A path to a file containing the expression ratios.
#' @param sep The separator between columns.
#' @param header Does the first data line contain column names?
#' @param ... Extra arguments passed on to [data.table::fread()].
#'
#' @return A data frame of expression values. Each row is for a locus. The locus
#'   identity is indicated in the first column. Remaining columns are samples.
#' @md
#' @export
read_ae_ratios <- function(file, sep = ",", header = TRUE, ...) {

  # Import the expression ratios (AE ratios), one row per locus.
  # First column is expected to be the SNP identifier, remaining columns are
  # the samples.
  dt <- data.table::fread(input = file, sep = sep, header = header, ...)

  # Fix the name of the first column
  data.table::setnames(dt, 1L, 'dae_snp')

  # Create key on `dae_snp`
  data.table::setkeyv(dt, 'dae_snp')

  dt
}
