% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_snp_pairs.R
\name{read_snp_pairs}
\alias{read_snp_pairs}
\title{Import pairs of DAE and DAEQTL candidate SNPs}
\usage{
read_snp_pairs(file, ...)
}
\arguments{
\item{file}{A path to a file.}

\item{...}{Extra arguments to be passed on to \code{\link[data.table:fread]{data.table::fread()}}.}
}
\value{
A data frame.
}
\description{
This function reads in a table of pairs of SNPs from a CSV file: the DAE SNP
and candidate DAEQTL SNP. One pair per row. The file must have a header with
the columns:
\itemize{
\item \code{dae_snp}: the DAE SNP.
\item \code{candidate_snp}: the candidate DAEQTL SNP.
\item \code{chromosome}: the chromosome where the pair is located.
\item \code{dae_snp_position}: the DAE SNP position on the chromosome.
\item \code{candidate_snp_position}: the candidate DAEQTL SNP position on the chromosome.
}
}
\seealso{
\code{\link[=create_snp_pairs]{create_snp_pairs()}}
}
