#' @keywords internal
daeqtl_mapping_ <- function(snp_pairs, zygosity, ae, fn = daeqtl_test, ...) {

  # Ensure that the data tables are all keyed.
  # `snp_pairs` keys: first two columns
  data.table::setkeyv(snp_pairs, colnames(snp_pairs)[1:2])
  # `zygosity` keys: first column
  data.table::setkeyv(zygosity, colnames(zygosity)[1])
  # `ae` keys: first column
  data.table::setkeyv(ae, colnames(ae)[1])

  n <- nrow(snp_pairs)
  for(i in seq_len(n)) {

    candidate_snp <- snp_pairs[[i, 'candidate_snp']]
    dae_snp <- snp_pairs[[i, 'dae_snp']]

    if(dae_snp %is_not_in% zygosity) {
      data.table::set(snp_pairs, i = i, j = 'pvalue', value = NA)
      data.table::set(snp_pairs, i = i, j = 'case', value = 0L)
      next
    }

    if(dae_snp %is_not_in% ae) {
      data.table::set(snp_pairs, i = i, j = 'pvalue', value = NA)
      data.table::set(snp_pairs, i = i, j = 'case', value = 0L)
      next
    }

    # `csnp_hom`: lgl indicating homozygous samples for the candidate snp
    csnp_hom <- is_hom(dt = zygosity, snp = candidate_snp)

    # `csnp_het`: lgl indicating heterozygous samples for the candidate snp
    csnp_het <- is_het(dt = zygosity, snp = candidate_snp)

    # `dsnp_het`: lgl indicating heterozygous samples for the dae snp
    dsnp_het <- is_het(dt = zygosity, snp = dae_snp)

    # `ae_hom`: dbl with allelic ratios of homozygous samples (candidate snp)
    ae_hom <- as.numeric(ae[dae_snp][, 1 := NULL][, dsnp_het & csnp_hom, with = FALSE])

    # `ae_het`: dbl with allelic ratios of heterozygous samples (candidate snp)
    ae_het <- as.numeric(ae[dae_snp][, 1 := NULL][, dsnp_het & csnp_het, with = FALSE])

    df <- fn(ae_hom = ae_hom[!is.na(ae_hom)], ae_het = ae_het[!is.na(ae_het)], ...)
    for (col in names(df)) data.table::set(snp_pairs, i = i, j = col, value = df[[col]])
  }

  return(snp_pairs[])
}

#' DAEQTL mapping
#'
#' @description
#' TODO
#'
#' @param snp_pairs A data frame with DAE SNP/candidate SNP pairs.
#' @param zygosity A data frame with zygosity levels of samples. First column
#'   should be the SNP identity, and remaining columns should refer to
#'   biological samples. These should be factors, and their levels should be
#'   `"hom"` (homozygous) and `"het"` (heterozygous).
#' @param ae A data frame with allelic expression (AE) ratios. Each row is for a
#'   locus. The locus identity is indicated in the first column. Remaining
#'   columns are samples.
#' @param fn A function implementing the statistical association approach.
#'   This function needs to have two named arguments `ae_hom` and `ae_het`.
#'   Extra arguments will be read from `...`.
#' @param ... Extra arguments passed on to the call of `test_fn()`.
#'
#' @return An updated version of the data frame `snp_pairs`. The update includes
#'   extra columns, typically `pvalue` and `case`.
#'
#' @md
#' @export
daeqtl_mapping <- daeqtl_mapping_
