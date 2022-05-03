#' Extract a vector of allelic expression ratios
#'
#' @description
#' `ae()` provides a handy shortcut to extract a vector of allelic
#' expression ratios for a DAE SNP on a set of samples.
#'
#' `ae_hom()` and `ae_het()` also extract a vector of allelic expression ratios
#' but for the corresponding homozygous or heterozygous samples of the candidate
#' SNP, respectively.
#'
#' @param dae_snp A string indicating the identifier of the DAE SNP.
#' @param ae A data frame of expression values. Each row is for a locus. The
#'   locus identity is indicated in the first column. Remaining columns are
#'   samples.
#' @param samples Either a logical or an integer vector indicating which samples
#' in the `ae` data table are to be selected.
#' @param drop_na Whether to drop `NA` values in the returned value.
#' @param candidate_snp A string indicating the identifier of the candidate SNP.
#' @param zygosity A data frame of zygosity levels: `"hom"` for homozygous or
#'   `"het"` for heterozygous. Each row is for a locus. The locus identity is
#'   indicated in the first column and named \code{snp}. Remaining columns are
#'   samples.
#'
#' @return A numeric vector of allelic ratios.
#'
#' @md
#'
#' @examples
#' # Let us load some dummy data
#' zygosity <- read_snp_zygosity(file = daeqtlr_example("zygosity.csv"))
#' ae <- read_ae_ratios(file = daeqtlr_example("ae.csv"))
#'
#' # Select all allelic expression ratios of rsX019
#' ae('rsX019', ae = ae, drop_na = FALSE)
#'
#' # Select only the first 5 samples
#' ae('rsX019', ae = ae, samples = 1:5, drop_na = FALSE)
#'
#' # Use a logical vector to select samples that meet a requirement, e.g.
#' # heterozygous samples only. Note that `is_het()` is useful here.
#' (het_samples_for_rsX019 <- is_het(zygosity = zygosity, snp = 'rsX019'))
#' ae('rsX019', ae = ae, samples = het_samples_for_rsX019, drop_na = FALSE)
#'
#' # If you want the allelic ratios for samples that are simultaneously:
#' # - heterozygous for rsX019
#' # - homozygous for rsX002
#' (hom_samples_for_rsX002 <- is_hom(zygosity = zygosity, snp = 'rsX002'))
#' ae(
#'   dae_snp = 'rsX019',
#'   ae = ae,
#'   samples = het_samples_for_rsX019 & hom_samples_for_rsX002,
#'   drop_na = FALSE
#'   )
#'
#' # Or more simply:
#' ae_hom('rsX019', 'rsX002', zygosity, ae, drop_na = FALSE)
#'
#'
#' @export
ae <- function(dae_snp, ae, samples = rep(TRUE, ncol(ae) - 1), drop_na = TRUE) {

  ae_vec <- as.numeric(ae[dae_snp][, 1 := NULL][, samples, with = FALSE])

  if (drop_na) {
    return(ae_vec[!is.na(ae_vec)])
  } else {
    return(ae_vec)
  }

}

#' @rdname ae
#' @export
ae_hom <- function(dae_snp, candidate_snp, zygosity, ae, drop_na = TRUE) {

  # `csnp_hom`: lgl indicating homozygous samples for the candidate snp
  csnp_hom <- is_hom(zygosity = zygosity, snp = candidate_snp)

  # `csnp_het`: lgl indicating heterozygous samples for the candidate snp
  csnp_het <- is_het(zygosity = zygosity, snp = candidate_snp)

  # `dsnp_het`: lgl indicating heterozygous samples for the dae snp
  dsnp_het <- is_het(zygosity = zygosity, snp = dae_snp)

  # `ae_hom`: dbl with allelic ratios of homozygous samples (candidate snp)
  # ae_hom <- as.numeric(ae[dae_snp][, 1 := NULL][, dsnp_het & csnp_hom, with = FALSE])
  ae_hom <- ae(dae_snp = dae_snp, ae = ae, samples = dsnp_het & csnp_hom, drop_na = drop_na)

  ae_hom
}

#' @rdname ae
#' @export
ae_het <- function(dae_snp, candidate_snp, zygosity, ae, drop_na = TRUE) {

  # `csnp_hom`: lgl indicating homozygous samples for the candidate snp
  csnp_hom <- is_hom(zygosity = zygosity, snp = candidate_snp)

  # `csnp_het`: lgl indicating heterozygous samples for the candidate snp
  csnp_het <- is_het(zygosity = zygosity, snp = candidate_snp)

  # `dsnp_het`: lgl indicating heterozygous samples for the dae snp
  dsnp_het <- is_het(zygosity = zygosity, snp = dae_snp)

  # `ae_het`: dbl with allelic ratios of heterozygous samples (candidate snp)
  # ae_het <- as.numeric(ae[dae_snp][, 1 := NULL][, dsnp_het & csnp_het, with = FALSE])
  ae_het <- ae(dae_snp = dae_snp, ae = ae, samples = dsnp_het & csnp_het)

  ae_het
}


