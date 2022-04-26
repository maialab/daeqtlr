#' Import pairs of DAE and DAEQTL candidate SNPs
#'
#' @description
#' This function reads in a table of pairs of SNPs from a CSV file: the DAE SNP
#' and candidate DAEQTL SNP. One pair per row. The file must have a header with
#' the columns:
#' - `dae_snp`: the DAE SNP.
#' - `candidate_snp`: the candidate DAEQTL SNP.
#' - `chromosome`: the chromosome where the pair is located.
#' - `dae_snp_position`: the DAE SNP position on the chromosome.
#' - `candidate_snp_position`: the candidate DAEQTL SNP position on the chromosome.
#'
#' @param file A path to a file.
#' @param ... Extra arguments to be passed on to [data.table::fread()].
#'
#' @return A data frame.
#'
#' @seealso [create_snp_pairs()]
#' @md
#' @export
read_snp_pairs <- function(file, ...) {

  col_names <-
    c(
      'dae_snp',
      'candidate_snp',
      'chromosome',
      'dae_snp_position',
      'candidate_snp_position'
    )

  col_types <-
    c('character', 'character', 'character', 'integer', 'integer')

  snp_pairs <-
    data.table::fread(
      input = file,
      colClasses = col_types,
      col.names = col_names,
      ...
    )

  snp_pairs
}
