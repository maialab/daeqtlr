#' Import SNP genomic positions
#'
#' This function reads in genomic positions from a text file. The data is
#' expected to be in tabular format. First column is expected to be the locus
#' identity (e.g. SNP identifier), second column is the human chromosome (e.g.,
#' 22 or X), and the third column is the position.
#'
#' @param file A path to a file containing the SNP genomic positions.
#' @param sep The separator between columns.
#' @param header Does the first data line contain column names?
#' @param ... Extra arguments passed on to [data.table::fread()].
#'
#' @return A data frame of genomic positions. Each row is for a locus. The locus
#'   identity is indicated in the first column (\code{snp}). Second column is
#'   the chromosome, and the third column is the position.
#' @md
#' @export
read_snp_gen_positions <- function(file, sep = ',', header = TRUE, ...) {

  genomic_positions <-
    data.table::fread(
      input = file,
      sep = sep,
      header = header,
      colClasses = c('character', 'character', 'integer'),
      col.names = c('snp', 'chromosome', 'position'),
      ...
    )

  genomic_positions
}
