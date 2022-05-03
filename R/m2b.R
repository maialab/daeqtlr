#' Convert M-values to Beta-values and back again
#'
#' - `m2b()` converts M-values to Beta-values.
#' - `b2m()` converts Beta-values to M-values.
#'
#' @param m A numeric vector of M-values.
#' @param b A numeric vector of Beta-values.
#' @param base The base with respect to which logarithms are computed. Defaults
#'   to `base = 2`.
#'
#' @md
#' @references Du, P. et al. _Comparison of Beta-value and M-value methods for
#'   quantifying methylation levels by microarray analysis_. BMC Bioinformatics
#'   11, (2010). \doi{10.1186/1471-2105-11-587}.
#'
#' @examples
#' # M-values and Beta-values have different domains
#' m2b(c(-Inf, 0, Inf))
#' b2m(c(0, 0.5, 1))
#'
#' # `m2b()` and `b2m()` are the inverse of one another
#' b2m(m2b(c(-Inf, 0, Inf)))
#'
#' m2b(b2m(c(0, 0.5, 1)))
#'
#' @export
m2b <- function(m, base = 2) {

  # base ^ m / (base ^ m + 1)

  # This rewrite works with extreme M-values, i.e. with both -Inf and Inf,
  # returning appropriately 0 and 1.
  1 / (1 + base ^ (-m))

}

#' @rdname m2b
#' @export
b2m <- function(b, base = 2) {

  logb(b / (1 - b), base = base)

}
