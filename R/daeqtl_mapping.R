#' @keywords internal
daeqtl_mapping_ <- function(snp_pairs, zygosity, ae) {

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

    if(dae_snp %is_not_in% snp_pairs) {
      datatable::set(snp_pairs, i = i, j = 'pvalue', value = NA)
      datatable::set(snp_pairs, i = i, j = 'case', value = 0L)
      break
    }

    if(dae_snp %is_not_in% ae) {
      datatable::set(snp_pairs, i = i, j = 'pvalue', value = NA)
      datatable::set(snp_pairs, i = i, j = 'case', value = 0L)
      break
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


    df <- daeqtl_test(ae_hom = ae_hom, ae_het = ae_het)
    for (col in names(df)) data.table::set(snp_pairs, i = i, j = col, value = df[[col]])
  }
}
