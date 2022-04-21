ae_lim <- function(ae) {

  ceiling(max(max(ae), 1))

}

#' @export
daeqtl_plot <- function(ae_hom,
                        ae_het,
                        dae_threshold = log2(1.5),
                        min_n_hom = 2L,
                        min_n_het = 2L,
                        dae_threshold_linetype = 'dashed',
                        xlab = 'Candidate SNP zygosity',
                        ylab = 'DAE SNP allelic expression (log-ratio)') {

  n_hom <- length(ae_hom)
  n_het <- length(ae_het)

  zygosity <- factor(c(rep('hom', n_hom), rep('het', n_het)), levels = c('hom', 'het'))
  ae <- c(ae_hom, ae_het)

  tbl <-
    tibble::tibble(zygosity = zygosity,
                   ae = ae)

  y_lim <- ae_lim(ae)

  ggplot2::ggplot(data = tbl, mapping = ggplot2::aes(x = zygosity, y = ae)) +
    ggplot2::scale_x_discrete(drop = FALSE) + # ensures all levels are kept
    ggplot2::geom_point() +
    ggplot2::ylim(-y_lim, y_lim) +
    ggplot2::geom_hline(yintercept = c(-dae_threshold, dae_threshold), linetype = dae_threshold_linetype) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)

}
