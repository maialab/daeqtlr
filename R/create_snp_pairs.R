# #' Create a table of SNP pairs
# #'
# #' This function creates a data table of SNP pairs to be used in DAEQTL mapping.
# #' Essentially, this function looks for neighbouring SNPs, i.e. within a
# #' genomic window (specified with `window_size`) and exports a file with those
# #' SNP pairs, in long format. Check the `vignette('snp-pairs')` to understand
# #' when to use this function.
# #'
# #' @param file A path to a file where the newly created data table of SNP pairs
# #' is to be exported.
# #' @param snp_table A data table of SNP to be used in DAEQTL mapping. Typically,
# #' this object is created with [create_snp_table()].
# #' @param window_size A genomic window size, in base pairs. Default is 500 kb.
# #' @param verbose Whether to be chatty about this function's progress.
# #'
# #' @return This function is run for its side effect: creating the file indicated
# #' in `file`.
# #'
# #' @md
# #' @importFrom rlang .data
# #' @export
# create_snp_pairs <-
#   function(file,
#            snp_table,
#            window_size = 500000L,
#            verbose = TRUE) {
#
#     dae_snp_tbl <-
#       dplyr::filter(snp_table, .data$is_dae_snp)
#
#     chromosomes <- unique(dae_snp_tbl$chromosome)
#     chromosomes <- chromosomes[!is.na(chromosomes)]
#
#     append <- FALSE
#
#     for (i in seq_along(chromosomes)) {
#       chr <- chromosomes[i]
#
#       # DAE SNPs at chromosome `chr`:
#       dae_snps_by_chr <-
#         dplyr::filter(dae_snp_tbl, .data$chromosome == chr)
#
#       # Candidate SNPs at chromosome `chr`:
#       candidate_snps_by_chr <-
#         dplyr::filter(snp_table, .data$chromosome == chr, .data$is_candidate_snp) %>%
#         dplyr::select(-c('is_dae_snp', 'is_candidate_snp'))
#
#       if (verbose) {
#         cli::cli_alert_info("Chromosome {chr}")
#         cli::cli_progress_bar("Writing", total = nrow(dae_snps_by_chr))
#       }
#       for (j in seq_along(dae_snps_by_chr$snp)) {
#         dae_snp <- dae_snps_by_chr[['snp']][j]
#         dae_snp_pos <- dae_snps_by_chr[['position']][j]
#
#         # Collect neighboring candidate SNPs around the DAE SNP
#         candidate_snps_by_chr_by_dae_snp <-
#           candidate_snps_by_chr %>%
#           dplyr::filter(.data$snp != dae_snp) %>% # Exclude self
#           dplyr::filter(abs(dae_snp_pos - .data$position) <= window_size) %>%
#           dplyr::arrange(.data$position) %>%
#           dplyr::rename(candidate_snp = .data$snp,
#                         candidate_snp_position = .data$position) %>%
#           dplyr::mutate(dae_snp = dae_snp,
#                         dae_snp_position = dae_snp_pos) %>%
#           dplyr::relocate(
#             .data$dae_snp,
#             .data$candidate_snp,
#             .data$chromosome,
#             .data$dae_snp_position,
#             .data$candidate_snp_position
#           )
#
#         data.table::fwrite(candidate_snps_by_chr_by_dae_snp,
#                            file = file,
#                            append = append)
#         append <- TRUE
#         if (verbose)
#           cli::cli_progress_update()
#       }
#     }
#     if (verbose)
#       cli::cli_alert_info("{file} written!")
#
#     return(file)
#   }

#' Create a table of SNP pairs
#'
#' This function creates a data table of SNP pairs to be used in DAEQTL mapping.
#' Essentially, this function looks for neighboring SNPs, i.e. within a genomic
#' window (specified with `window_size`) and creates a data table with those SNP
#' pairs, in long format. Check the `vignette('snp-pairs')` to understand when
#' you might need to use this function.
#'
#' @param snp_table A data table of SNP to be used in DAEQTL mapping. Typically,
#'   this object is created with [create_snp_table()]. See
#'   `vignette('snp-pairs')` for context.
#' @param window_size A genomic window size, in base pairs. Default is 500 kb,
#'   up and down the DAE SNP, i.e. a window of 1 Mb around the DAE SNP.
#'
#' @return A data table.
#'
#' @md
#' @importFrom rlang .data
#' @export
create_snp_pairs <- function(snp_table, window_size = 500000L) {

  # Because of NSE: https://github.com/Rdatatable/data.table/issues/850#issuecomment-259466153
  position <- NULL
  is_dae_snp <- is_candidate_snp <- NULL
  dae_snp <- candidate_snp <- NULL

  dsnp_dt <- snp_table[(is_dae_snp)][, c("is_dae_snp", "is_candidate_snp") := NULL]
  data.table::setnames(dsnp_dt, 'snp', 'dae_snp')
  dsnp_dt[, c('start', 'end') := list(position - window_size, position + window_size)]
  data.table::setnames(dsnp_dt, 'position', 'dae_snp_position')

  csnp_dt <- snp_table[(is_candidate_snp)][, c("is_dae_snp", "is_candidate_snp") := NULL]
  data.table::setnames(csnp_dt, 'snp', 'candidate_snp')
  csnp_dt[, c('start', 'end') := list(position, position)]
  data.table::setnames(csnp_dt, 'position', 'candidate_snp_position')

  data.table::setkeyv(dsnp_dt, c('chromosome', 'start', 'end'))
  data.table::setkeyv(csnp_dt, c('chromosome', 'start', 'end'))

  dt <- data.table::foverlaps(csnp_dt, dsnp_dt, type = 'within', nomatch = NULL)
  dt2 <- dt[, c('dae_snp', 'candidate_snp', 'chromosome', 'dae_snp_position', 'candidate_snp_position')]
  data.table::setkeyv(dt2, c('chromosome', 'dae_snp_position'))

  dt2
}
