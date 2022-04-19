`%is_in%` <- function(snp, dt) {
  # !(dt[.(snp), mult = "first", nomatch = 0L, .N] == 0)
  !(dt[snp, mult = "first", nomatch = NULL, .N] == 0L)
}
