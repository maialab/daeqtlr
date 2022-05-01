
<!-- README.md is generated from README.Rmd. Please edit that file -->

# daeqtlr

<!-- badges: start -->
<!-- badges: end -->

The goal of `{daeqtlr}` is to perform differential allelic expression
trait loci (DAEQTL) mapping.

Differential allelic expression (DAE) can be treated as a quantitative
trait. Marker genotypes collected from the same set of individuals can,
in turn, be used in testing statistical associations with DAE. This
package provides routines for performing such a workflow, i.e. DAEQTL
mapping.

## Installation

You can install the development version of `{daeqtlr}` like so:

``` r
# install.packages("remotes")
remotes::install_github("maialab/daeqtlr")
```

## Usage

Here’s a quick example on how to perform DAEQTL mapping with
`{daeqtlr}`. For more details about this example read
`vignette("daeqtlr")`.

``` r
library(daeqtlr)

# SNP pairs to be tested for statistical association
snp_pairs <- read_snp_pairs(file = daeqtlr_example("snp_pairs.csv"))

# Zygosity levels (homozygous or heterozygous) for genotyped SNPs
zygosity <- read_snp_zygosity(file = daeqtlr_example("zygosity.csv"))

# Allelic expression (AE) ratios of DAE SNPs
ae <- read_ae_ratios(file = daeqtlr_example("ae.csv"))

# `snp_pairs` is updated by reference
daeqtl_mapping(snp_pairs = snp_pairs, zygosity = zygosity, ae = ae)

# Here's the first 30 SNP pairs tested
snp_pairs[1:30, -c('dae_snp_position', 'candidate_snp_position')]
#>     dae_snp candidate_snp chromosome    pvalue case
#>  1:  rsX019        rsX002         19 0.8331668    4
#>  2:  rsX019        rsX003         19 0.1348031    4
#>  3:  rsX019        rsX005         19 0.9280741    4
#>  4:  rsX019        rsX010         19 0.9280741    4
#>  5:  rsX019        rsX014         19 0.1348031    4
#>  6:  rsX019        rsX016         19        NA    2
#>  7:  rsX019        rsX017         19        NA    2
#>  8:  rsX019        rsX018         19 0.9255817    4
#>  9:  rsX019        rsX022         19 0.3994533    4
#> 10:  rsX019        rsX023         19 0.3736310    4
#> 11:  rsX019        rsX024         19 0.8813157    4
#> 12:  rsX019        rsX025         19        NA    1
#> 13:  rsX019        rsX026         19        NA    1
#> 14:  rsX019        rsX043         19        NA    1
#> 15:  rsX019        rsX045         19        NA    2
#> 16:  rsX019        rsX046         19        NA    1
#> 17:  rsX019        rsX047         19        NA    2
#> 18:  rsX019        rsX048         19        NA    2
#> 19:  rsX019        rsX050         19        NA    2
#> 20:  rsX019        rsX051         19        NA    2
#> 21:  rsX019        rsX053         19        NA    2
#> 22:  rsX019        rsX054         19        NA    2
#> 23:  rsX019        rsX060         19 0.1348031    4
#> 24:  rsX019        rsX064         19        NA    1
#> 25:  rsX020        rsX008          1 0.1226391    4
#> 26:  rsX020        rsX011          1 0.1003482    3
#> 27:  rsX020        rsX012          1 0.1003482    3
#> 28:  rsX020        rsX013          1        NA    1
#> 29:  rsX020        rsX015          1 0.3492677    4
#> 30:  rsX020        rsX021          1 0.1814492    3
#>     dae_snp candidate_snp chromosome    pvalue case
```
