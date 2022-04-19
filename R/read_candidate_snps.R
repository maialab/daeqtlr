#' Import DAEQTL candidates
#'
#' This function reads in a list of candidate SNPs to being DAE quantitative
#' trait loci, from a text file. One SNP per row. By default, it expects a
#' header, which will be ignored. If the first line is immediately a SNP,
#' then use \code{header = FALSE} to not ignore it.
#'
#' @param file A path to a file.
#' @param header Is the first line a line with a column name?
#' @param distinct Remove any duplicate SNPs?
#'
#' @return A one-column data frame. This single column is named \code{snp}.
#'   Rows are the candidate SNPs read from the text file.
#'
#' @export
read_candidate_snps <- function(file, header = TRUE, distinct = TRUE) {
  read_snps(file, colname = 'candidate_snp', header = header, distinct = distinct)
}
