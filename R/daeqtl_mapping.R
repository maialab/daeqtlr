#' @keywords internal
daeqtl_mapping_ <- function(snp_pairs, zygosity, ae, fn = daeqtl_test, ..., .extra_cols = 2L) {

  # Ensure that the data tables are all keyed.
  # `snp_pairs` keys: first two columns
  data.table::setkeyv(snp_pairs, colnames(snp_pairs)[1:2])
  # `zygosity` keys: first column
  data.table::setkeyv(zygosity, colnames(zygosity)[1])
  # `ae` keys: first column
  data.table::setkeyv(ae, colnames(ae)[1])

  n <- nrow(snp_pairs)
  # `setalloccol` is needed because of `future.apply::future_lapply()`,
  # otherwise https://github.com/Rdatatable/data.table/issues/5376.
  data.table::setalloccol(snp_pairs, .extra_cols)

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

    ae_hom <- ae_hom(dae_snp, candidate_snp, zygosity, ae)
    ae_het <- ae_het(dae_snp, candidate_snp, zygosity, ae)

    df <- fn(ae_hom = ae_hom, ae_het = ae_het, ...)

    for (col in names(df)) data.table::set(snp_pairs, i = i, j = col, value = df[[col]])
  }

  return(invisible(snp_pairs[]))
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
#' @param .extra_cols The number of extra columns you are creating with `fn`.
#' @param .n_cores The number of chunks to divide `snp_pairs` for parallel
#'   processing.
#'
#' @return An updated version of the data frame `snp_pairs`. The update includes
#'   extra columns, typically `pvalue` and `case`.
#'
#' @md
#' @export
daeqtl_mapping <-
  function(snp_pairs,
           zygosity,
           ae,
           fn = daeqtl_test,
           ...,
           .extra_cols = 2L,
           .n_cores = 1L) {

    snp_pairs_lst <-
      split(snp_pairs, split_index(nrow(snp_pairs), .n_cores))

    for (i in seq_along(snp_pairs_lst)) {
      data.table::setkeyv(snp_pairs_lst[[i]], c('dae_snp', 'candidate_snp'))
    }

    future::plan(future::multisession)

    res <- future.apply::future_lapply(
      snp_pairs_lst,
      FUN = daeqtl_mapping_,
      zygosity = zygosity,
      ae = ae,
      fn = fn,
      ...,
      .extra_cols = .extra_cols
    )

    future::plan(future::sequential)

    # Put together the data table from the list of data tables
    mapping_dt <- data.table::rbindlist(res)

    data.table::setkeyv(mapping_dt, cols = c('dae_snp', 'candidate_snp'))
    return(mapping_dt)
}

#' @export
daeqtl_mapping2 <-
  function(snp_pairs,
           zygosity,
           ae,
           fn = daeqtl_test,
           ...,
           .extra_cols = 2L,
           .n_cores = 1L) {

    snp_pairs_lst <-
      split(snp_pairs, split_index(nrow(snp_pairs), .n_cores))

    for (i in seq_along(snp_pairs_lst)) {
      data.table::setkeyv(snp_pairs_lst[[i]], c('dae_snp', 'candidate_snp'))
    }

    res <- future.apply::future_lapply(
      snp_pairs_lst,
      FUN = daeqtl_mapping_,
      zygosity = zygosity,
      ae = ae,
      fn = fn,
      ...,
      .extra_cols = .extra_cols
    )

    # Put together the data table from the list of data tables
    mapping_dt <- data.table::rbindlist(res)

    data.table::setkeyv(mapping_dt, cols = c('dae_snp', 'candidate_snp'))
    return(mapping_dt)
  }
