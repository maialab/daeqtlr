#' Create a table of SNP pairs
#'
#' This function creates a data table of SNP pairs to be used in DAEQTL mapping.
#' Essentially, this function looks for neighbouring SNPs, i.e. within a
#' genomic window (specified with `window_size`) and exports a file with those
#' SNP pairs, in long format.
#'
#' @param file A path to a file where the newly created data table of SNP pairs
#' is to be exported.
#' @param snp_table A data table of SNP to be used in DAEQTL mapping. Typically,
#' this object is created with [create_snp_table()].
#' @param window_size A genomic window size, in base pairs. Default is 500 kb.
#' @param verbose Whether to be chatty about this function's progress.
#'
#' @return This function is run for its side effect: creating the file indicated
#' in `file`.
#'
#' @md
#' @importFrom rlang .data
#' @export
create_snp_pairs <-
  function(file,
           snp_table,
           window_size = 500000L,
           verbose = TRUE) {

    dae_snp_tbl <-
      dplyr::filter(snp_table, .data$is_dae_snp)

    chromosomes <- unique(dae_snp_tbl$chromosome)
    chromosomes <- chromosomes[!is.na(chromosomes)]

    append <- FALSE

    for (i in seq_along(chromosomes)) {
      chr <- chromosomes[i]

      # DAE SNPs at chromosome `chr`:
      dae_snps_by_chr <-
        dplyr::filter(dae_snp_tbl, .data$chromosome == chr)

      # Candidate SNPs at chromosome `chr`:
      candidate_snps_by_chr <-
        dplyr::filter(snp_table, .data$chromosome == chr, .data$is_candidate_snp) %>%
        dplyr::select(-c('is_dae_snp', 'is_candidate_snp'))

      if (verbose) {
        cli::cli_alert_info("Chromosome {chr}")
        cli::cli_progress_bar("Writing", total = nrow(dae_snps_by_chr))
      }
      for (j in seq_along(dae_snps_by_chr$snp)) {
        dae_snp <- dae_snps_by_chr[['snp']][j]
        dae_snp_pos <- dae_snps_by_chr[['position']][j]

        # Collect neighboring candidate SNPs around the DAE SNP
        candidate_snps_by_chr_by_dae_snp <-
          candidate_snps_by_chr %>%
          dplyr::filter(.data$snp != dae_snp) %>% # Exclude self
          dplyr::filter(abs(dae_snp_pos - .data$position) <= window_size) %>%
          dplyr::arrange(.data$position) %>%
          dplyr::rename(candidate_snp = .data$snp,
                        candidate_snp_position = .data$position) %>%
          dplyr::mutate(dae_snp = dae_snp,
                        dae_snp_position = dae_snp_pos) %>%
          dplyr::relocate(
            .data$dae_snp,
            .data$candidate_snp,
            .data$chromosome,
            .data$dae_snp_position,
            .data$candidate_snp_position
          )

        data.table::fwrite(candidate_snps_by_chr_by_dae_snp,
                           file = file,
                           append = append)
        append <- TRUE
        if (verbose)
          cli::cli_progress_update()
      }
    }
    if (verbose)
      cli::cli_alert_info("{file} written!")

    return(file)
  }
