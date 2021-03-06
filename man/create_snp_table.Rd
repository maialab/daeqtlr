% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_snp_table.R
\name{create_snp_table}
\alias{create_snp_table}
\title{Create a SNP table}
\usage{
create_snp_table(snp_gen_pos, dae_snps, candidate_snps)
}
\arguments{
\item{snp_gen_pos}{A data frame of three columns: \code{snp},
\code{chromosome} and \code{position}.}

\item{dae_snps}{A data frame of one column containing the DAE SNPs.}

\item{candidate_snps}{A data frame of one column containing the DAEQTL candidate SNPs.}
}
\value{
This function adds two columns to \code{snp_gen_pos} in-place, so the
object passed in \code{snp_gen_pos} will be modified after the call to this
function.
}
\description{
This function assembles a data frame of all SNPs involved in the DAEQTL
mapping. Besides the genomic position for each SNP, it annotates each SNP as
a DAE SNP or a DAEQTL candidate SNP: columns \code{is_dae_snp} and
\code{is_candidate_snp}.

The returned data frame is useful as an intermediate data structure of the
workflow, namely when looking for candidate neighbouring SNPs. Check the
\code{vignette('snp-pairs')} to understand when to use this function.
}
