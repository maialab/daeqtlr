#' Convert M-values to Beta-values
#'
#' This function converts M-values to Beta-values.
#'
#' @param m A numeric vector of M-values.
#' @param base The base with respect to which logarithms are computed. Defaults
#'   to `base = 2`.
#'
#' @md
#' @examples
#' m2b(c(-Inf, 0, Inf))
#'
#' @seealso [daeqtlr::b2m()] for the reverse operation.
#' @export
m2b <- function(m, base = 2) {

  # base ^ m / (base ^ m + 1)

  # This rewrite works with extreme M-values, i.e. with both -Inf and Inf,
  # returning appropriately 0 and 1.
  1 / (1 + base ^ (-m))

}

#' Convert Beta-values to M-values
#'
#' This function converts Beta-values to M-values.
#'
#' @param b A numeric vector of Beta-values.
#' @param base The base with respect to which logarithms are computed. Defaults
#'   to `base = 2`.
#'
#' @md
#' @examples
#' b2m(c(0, 0.5, 1))
#'
#' @seealso [daeqtlr::m2b()] for the reverse operation.
#' @export
b2m <- function(b, base = 2) {

  logb(b / (1 - b), base = base)

}
