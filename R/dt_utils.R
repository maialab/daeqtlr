#' Who are the homozygous/heterozygous samples?
#'
#' This function determines who are the homozygous or heterozygous samples given
#' a data table of zygosity levels and a SNP of interest. Samples whose zygosity
#' level is `NA` will be returned as `FALSE`.
#'
#' @param zygosity A data frame of zygosity levels: `"hom"` for homozygous or
#'   `"het"` for heterozygous. Each row is for a locus. The locus identity is
#'   indicated in the first column and named \code{snp}. Remaining columns are
#'   samples.
#' @param snp String with SNP identifier.
#' @param na_as_false Whether to return `FALSE` when the zygosity level is `NA`.
#'
#' @return A logical vector.
#'
#' @md
#' @examples
#' # Let us start by reading in an example data set with zygosity levels
#' zygosity <- read_snp_zygosity(file = daeqtlr_example("zygosity.csv"))
#'
#' # Checking out SNP rsX005
#' zygosity['rsX005']
#' is_hom(zygosity, 'rsX005')
#' is_het(zygosity, 'rsX005')
#'
#' # Translate the logical vector to sample names
#' # Note that first column is excluded because it is the SNP identifier col
#' # Homozygous samples
#' (homs <- colnames(zygosity)[-1][is_hom(zygosity, 'rsX005')])
#'
#' # Heterozygous samples
#' (hets <- colnames(zygosity)[-1][is_het(zygosity, 'rsX005')])
#'
#' # Some samples are neither homozygous nor heterozygous because of NAs
#' # Note the `- 1` because the first column of `zygosity` is the SNP identifier.
#' (ncol(zygosity) - 1) - length(c(homs, hets))
#'
#' # The samples whose zygosity is NA are:
#' setdiff(colnames(zygosity)[-1], c(homs, hets))
#'
#' @export
is_hom <- function(zygosity, snp, na_as_false = TRUE) {
  x <- zygosity[snp, mult = "first", nomatch = 0L][, 1 := NULL]
  if (na_as_false) {
    lgl <- as.vector(x == 'hom' & !is.na(x))
  } else {
    lgl <- as.vector(x == 'hom')
  }
  return(lgl)
}

#' @rdname is_hom
#' @export
is_het <- function(zygosity, snp, na_as_false = TRUE) {
  x <- zygosity[snp, mult = "first", nomatch = 0L][, 1 := NULL]

  if (na_as_false) {
    lgl <- as.vector(x == 'het' & !is.na(x))
  } else {
    lgl <- as.vector(x == 'het')
  }
  return(lgl)
}
