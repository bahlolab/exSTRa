\name{exSTRa}
\alias{exstra}
\docType{package}
\title{
Expanded STR algorithm
}
\description{
Detecting expansions with paired-end Illumina sequencing data
}
\details{
\tabular{ll}{
Package: \tab exSTRa\cr
Type: \tab Package\cr
Version: \tab 0.9.1\cr
Date: \tab 2022-08-18\cr
License: \tab GPL-2\cr
}
For an overview of how to use the package, see TODO. 
}
\author{
Maintainer: \email{tankard.r@wehi.edu.au}
}
\references{
Detecting known repeat expansions with standard protocol next generation sequencing, towards developing a single screening test for neurological repeat expansion disorders.
Rick M. Tankard, Martin B. Delatycki, Paul J. Lockhart, Melanie Bahlo.
bioRxiv 157792; June 2017. doi: \url{https://doi.org/10.1101/157792}
}
\keyword{ package }
\seealso{
Data must first be prepared from BAM/CRAM files with the Perl module 
\url{https://github.com/bahlolab/Bio-STR-exSTRa}.
}
\examples{
exstra_wgs_pcr_2["HD"]

plot(exstra_wgs_pcr_2["HD"])

}
