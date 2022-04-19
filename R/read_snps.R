read_snps <- function(file, colname, header = TRUE, distinct = TRUE) {

  skip <- if(header) 1 else 0

  tbl <- readr::read_delim (
    file = file,
    delim = " ",
    col_names = colname,
    col_types = "c",
    skip = skip
  )

  if (distinct)
    tbl <- dplyr::distinct(tbl)

  tbl
}
