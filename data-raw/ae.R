# This code makes homozygous samples for the DAE SNP to be not available (NA)
# in the `ae` data set.

library(daeqtlr)
library(here)

# Read original allelic ratios: this file contains allelic ratios for the
# homozygous samples which is non-sensical, and only happens because these
# values are originally derived from a microarray experiment using
# llumina’s HumanExon510s-Duo genotyping BeadChip platform. These arrays probe
# the two alleles and give intensities for the two alleles regardless of
# zygosity. Therefore, for these homozygous samples these are non-sensical.
#
ae <- read_ae_ratios(file = here('data-raw', 'ae.csv'))

# Read zygosity levels
zygosity <- read_snp_zygosity(file = daeqtlr_example('zygosity.csv'))

# Make a matrix of zygosity levels, column `snp` dropped and each row
# corresponding to the same SNP in `ae`
zygosity_m <- zygosity[(ae$dae_snp)] %>%
  dplyr::select(-'snp') %>%
  as.matrix()

# Make a matrix out of the `ae` data table (drop `dae_snp` column)
ae_m <-
  ae %>%
  dplyr::select(-'dae_snp') %>%
  as.matrix()

# Set those allelic ratios that correspond to homozygous samples to NA
ae_m[zygosity_m == 'hom'] <- NA

# Generate new `ae` with updated values
ae <- data.frame(dae_snp = ae$dae_snp, ae_m)

# Finally export to inst/extdata
write.csv(
  ae,
  file = here('inst/extdata', 'ae.csv'),
  quote = FALSE,
  row.names = FALSE
)
