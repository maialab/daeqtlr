is_hom <- function(dt, snp) {
  x <- dt[snp, mult = "first", nomatch = 0L][, 1 := NULL]
  as.vector(x == 'hom' & !is.na(x))
}

is_het <- function(dt, snp) {
  x <- dt[snp, mult = "first", nomatch = 0L][, 1 := NULL]
  as.vector(x == 'het' & !is.na(x))
}
