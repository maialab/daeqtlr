#' Create a split index
#'
#' This function creates an integer vector that can be used as a factor for
#' grouping with `split()`.
#'
#' @param n Number of elements of the vector.
#' @param m Number of groups (or partitions).
#'
#' @return An integer vector.
#'
#' @keywords internal
split_index <- function(n, m) {

  if(n < 1) stop('`n` must be great than zero.')
  if(m < 1) stop('`m` must be great than zero.')
  if(m > n) stop('`m` must be less than or equal to `n`.')

  return(ceiling(seq_len(n) / (n / m)))

}
