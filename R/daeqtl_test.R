#' DAEQTL mapping approach
#'
#' @description
#' This function implements daeqtlr's default mapping approach, i.e. the mapping
#' of SNPs associated with DAE. The mapping approach here implemented takes into
#' consideration the pattern of the allelic expression (AE) ratios' distribution
#' displayed at each DAE SNP (depicted below), as this is dependent on the
#' linkage disequilibrium between the DAE SNP and the candidate SNP.
#'
#' \if{html}{\figure{mapping_approach.svg}{Mapping approach}}
#' \if{latex}{\figure{mapping_approach.png}{options: width=0.5in}}
#'
#' @param ae_hom Numeric vector of AE ratios of the DAE SNP. Each element of the
#'   vector refers to a sample that is homozygous (\code{hom}) for the candidate
#'   SNP.
#' @param ae_het Numeric vector of AE ratios of the DAE SNP. Each element of the
#'   vector refers to a sample that is heterozygous (\code{het}) for the
#'   candidate SNP.
#' @param dae_threshold An allelic expression (AE) threshold (in log-scale). A
#'   sample having showing an absolute AE greater than `dae_threshold` is
#'   considered "technically" \emph{differential allelic expressed}, meaning
#'   that the imbalance observed is not below the limit of detection. Adjustment
#'   made to this parameter should depend on the experimental assay sensitivity.
#' @param min_n_hom Minimum number of samples in the homozygous group to be
#'   considered elligible for statistical testing.
#' @param min_n_het Minimum number of samples in the heterozygous group to be
#'   considered elligible for statistical testing.
#'
#' @return A data frame of two columns:
#' \describe{
#' \item{pvalue}{The p-value associated with the statistical test, if performed; otherwise, \code{NA}.}
#' \item{case}{One of the four possible cases depicted in the figure above.
#' The identified case depends on the number of samples for each group and on
#' the pattern of the allelic expression (AE) ratios of the DAE SNP:
#' \describe{
#'   \item{`1`}{**Case 1**: the number of samples for which the candidate SNP is
#'   heterozygous and DAE SNP expression is available is below the minimum
#'   eligibility criterion \code{min_n_het}.}
#'   \item{`2`}{**Case 2**: (i) the number of samples for which the candidate SNP is
#'   homozygous is below \code{min_n_hom}, and also, (ii) the values of the DAE
#'   SNP AE ratios are not all either below or above the \code{dae_threshold}.}
#'   \item{`3`}{**Case 3**: the number of samples for which the candidate SNP is
#'   homozygous is below \code{min_n_hom}, however, the values of the DAE
#'   SNP AE ratios are all either below or above the \code{dae_threshold}. This
#'   case is followed by a one-sample Wilcox test whose null hypothesis is that
#'   the AE ratios are zero.}
#'   \item{`4`}{**Case 4**: when both a minimum number of heterozygous and
#'   homozygous samples for the candidate SNP are available, then a Wilcox test
#'   is applied that compares the two groups. The AE ratios are first
#'   transformed to absolute values. This is because we want to test departure
#'   from zero in either direction. The null hypothesis is that absolute AE
#'   ratios for the heterozygous group is less than or equal to the homozygous
#'   group.}
#' }
#' }
#' }
#'
#' @md
#' @export
daeqtlr_mapping <-
  function(ae_hom,
           ae_het,
           dae_threshold = log2(1.5),
           min_n_hom = 2L,
           min_n_het = 2L
  ) {

    n_hom <- length(ae_hom)
    n_het <- length(ae_het)

    # Case 1
    if(n_het < min_n_het) {
      return(data.frame(pvalue = NA_real_, case = 1L))
    }

    if(n_hom < min_n_hom) {
      # Case 2
      if (!(all(ae_het >= dae_threshold) ||
            all(ae_het <= -dae_threshold))) {
        return(data.frame(pvalue = NA_real_, case = 2L))
      } else { # Case 3
        wc_test <- stats::wilcox.test(x = ae_het, alternative = 'two.sided', mu = 0, exact = FALSE)
        return(data.frame(pvalue = wc_test$p.value, case = 3L))
      }
    }

    # Case 4
    wc_test <- stats::wilcox.test(x = abs(ae_het), y = abs(ae_hom), alternative = 'greater', exact = FALSE)
    return(data.frame(pvalue = wc_test$p.value, case = 4L))
  }
