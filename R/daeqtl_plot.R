ae_lim <- function(ae) {

  ae <- ae[is.finite(ae)]
  ceiling(max(max(abs(ae)), 1))

}

#' DAEQTL plot
#'
#' A simple wrapper around a ggplot2 plot showing allelic expression ratios for
#' two groups of samples: homozygous and heterozygous. This is a handy function
#' to quickly inspect the data underlying each DAEQTL mapping test.
#'
#' @param ae_hom Numeric vector of AE ratios of the DAE SNP. Each element of the
#'   vector refers to a sample that is homozygous (\code{hom}) for the candidate
#'   SNP.
#' @param ae_het Numeric vector of AE ratios of the DAE SNP. Each element of the
#'   vector refers to a sample that is heterozygous (\code{het}) for the
#'   candidate SNP.
#' @param dae_threshold An allelic expression (AE) threshold (in log-scale). A
#'   sample showing an absolute AE greater than `dae_threshold` is considered
#'   "technically" \emph{differential allelic expressed}, meaning that the
#'   imbalance observed is not below the limit of detection. Adjustment made to
#'   this parameter should depend on the experimental assay sensitivity.
#' @param min_n_hom Minimum number of samples in the homozygous group to be
#'   considered eligible for statistical testing.
#' @param min_n_het Minimum number of samples in the heterozygous group to be
#'   considered eligible for statistical testing.
#' @param dae_threshold_linetype Line type for the DAE threshold. See [Line
#'   type](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#lines) for
#'   options.
#' @param xlab Title of the x-axis.
#' @param ylab Title of the y-axis.
#' @param ylim A two element vector specifying the y-axis limits.
#'
#' @examples
#' ae_hom <- c(1.9, 2.1, 2 , 1.5, 1.4)
#' ae_het <- c(0.3, 0.6, 0.7)
#' daeqtl_plot(ae_hom = ae_hom, ae_het = ae_het)
#' @md
#' @importFrom rlang .data
#' @export
daeqtl_plot <- function(ae_hom,
                        ae_het,
                        dae_threshold = log2(1.5),
                        min_n_hom = 2L,
                        min_n_het = 2L,
                        dae_threshold_linetype = 'dashed',
                        xlab = 'Candidate SNP zygosity',
                        ylab = 'DAE SNP allelic expression (log-ratio)',
                        ylim = NULL) {

  n_hom <- length(ae_hom)
  n_het <- length(ae_het)

  zygosity <- factor(c(rep('hom', n_hom), rep('het', n_het)), levels = c('hom', 'het'))
  ae <- c(ae_hom, ae_het)

  tbl <-
    tibble::tibble(zygosity = zygosity,
                   ae = ae) %>%
    dplyr::filter(is.finite(.data$ae))

  ylim <- `if`(is.null(ylim), c(-ae_lim(ae), ae_lim(ae)), ylim)

  ggplot2::ggplot(data = tbl, mapping = ggplot2::aes(x = zygosity, y = ae)) +
    ggplot2::scale_x_discrete(drop = FALSE) + # ensures all levels are kept
    ggbeeswarm::geom_beeswarm(na.rm = TRUE) +
    ggplot2::ylim(ylim) +
    ggplot2::geom_hline(yintercept = c(-dae_threshold, dae_threshold), linetype = dae_threshold_linetype) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)

}
