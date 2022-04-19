#' Convert genotypes to zygosity
#'
#' @description
#' This function takes a data table with genotypes (`"AA"`, `"AB"`, `"BA"` or
#' `"BB"`) and converts them to their corresponding zygosity, i.e., to
#' homozygous (`"hom"`) or heterozygous (`"het"`).
#'
#' The argument `genotypes` is expected to be a data table whose first column is
#' the SNP identifier, and remaining columns refer to samples on which the SNP
#' has been genotyped. All columns except the first should be factors, and
#' should contain the levels `"AA"`, `"AB"`, `"BA"` and `"BB"`.
#'
#' Note that this function changes its input by reference, meaning that the
#' variable passed as `genotypes` will be modified after running `as_zygous()`.
#' This is for performance reasons.
#'
#' @param genotypes A data table whose first column is the SNP identifier (a
#'   character vector), and remaining columns refer to samples on which the SNP
#'   has been genotyped. All columns except the first should be factors, and
#'   should contain the levels `"AA"`, `"AB"`, `"BA"` and `"BB"`.
#'
#' @return A data table with the same structure as the data table passed in
#'   `genotypes`. The only difference is that the sample columns (all except the
#'   first) have their levels recoded to `"hom"` (homozygous) or `"het"`
#'   (heterozygous).
#'
#' @seealso [read_snp_genotypes()]
#' @md
#' @examples
#' levels <- c("AA", "AB", "BA", "BB")
#' (df <- data.frame(
#'   snp = c('rs123', 'rs456'),
#'   sample_01 = factor(c('AB', 'BA'), levels),
#'   sample02 = factor(c('AA', 'BB'), levels)
#' ))
#' as_zygous(df)
#' df
#'
#' @export
as_zygous <- function(genotypes) {

  # First column is the SNP identifier, remaining columns are samples.
  cols <- colnames(genotypes)[-1]

  # Recode genotypes to zygosity levels
  for (j in cols)
    data.table::set(
      genotypes,
      j = j,
      value = forcats::fct_recode(
        genotypes[[j]],
        hom = 'AA',
        hom = 'BB',
        het = 'AB',
        het = 'BA'
      )
    )

  genotypes
}
