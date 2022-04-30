#' Import SNP zygosity
#'
#' This function reads in SNP genotypes from a text file and converts them to
#' zygosity levels. The data is expected to be in tabular format. First column
#' is expected to be the locus identity (e.g. SNP identifier), and remaining
#' columns are the samples. The genotypes are expected to be encoded as `"AA"`,
#' `"AB"`, `"BA"` or `"BB"`.
#'
#' @param file A path to file.
#' @param sep The separator between columns.
#' @param header Does the first data line contain column names?
#' @param ... Extra arguments to be passed on to [data.table::fread()].
#'
#' @return A data frame of zygosity levels: `"hom"` for homozygous or `"het"`
#'   for heterozygous. Each row is for a locus. The locus identity is indicated
#'   in the first column and named \code{snp}. Remaining columns are samples.
#'
#' @md
#' @export
read_snp_zygosity <- function(file, sep = ',', header = TRUE, ...) {

  zygosity <- data.table::fread(input = file, sep = sep, header = header, ...)
  data.table::setnames(zygosity, 1L, "snp")
  data.table::setkeyv(zygosity, 'snp')

  # First column is the SNP identifier, remaining columns are samples.
  cols <- colnames(zygosity)[-1]

  # Convert genotype strings to factors
  for (j in cols)
    data.table::set(zygosity,
                    j = j,
                    value = factor(zygosity[[j]], levels = c('hom', 'het')))

  zygosity
}
