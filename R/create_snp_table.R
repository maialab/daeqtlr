#' Create a SNP table
#'
#' @description
#' This function assembles a data frame of all SNPs involved in the DAEQTL
#' mapping. Besides the genomic position for each SNP, it annotates each SNP as
#' a DAE SNP or a DAEQTL candidate SNP: columns \code{is_dae_snp} and
#' \code{is_candidate_snp}.
#'
#' The returned data frame is useful as an intermediate data structure of the
#' workflow, namely when looking for candidate neighbouring SNPs. Check the
#' `vignette('snp-pairs')` to understand when to use this function.
#'
#' @param snp_gen_pos A data frame of three columns: \code{snp},
#'   \code{chromosome} and \code{position}.
#' @param dae_snps A data frame of one column containing the DAE SNPs.
#' @param candidate_snps A data frame of one column containing the DAEQTL candidate SNPs.
#'
#' @return This function adds two columns to `snp_gen_pos` in-place, so the
#'   object passed in `snp_gen_pos` will be modified after the call to this
#'   function.
#'
#' @md
#' @export
create_snp_table <- function(snp_gen_pos, dae_snps, candidate_snps) {

  data.table::setkeyv(snp_gen_pos, 'snp')

  # boilerplate code for data.table
  is_dae_snp <- 'is_dae_snp'
  is_candidate_snp <- 'is_candidate_snp'

  snp_gen_pos[, (is_dae_snp) := get('snp') %in% dae_snps[[1]]]
  snp_gen_pos[, (is_candidate_snp) := get('snp') %in% candidate_snps[[1]]]

  # The `[]` is needed to force a print:
  # https://stackoverflow.com/questions/33195362/data-table-is-not-displayed-on-first-call-after-being-modified-in-a-function
  snp_gen_pos[]
}
