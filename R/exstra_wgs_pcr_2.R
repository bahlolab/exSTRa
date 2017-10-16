#' Repeat score data from cohort WGS_PCR_2.
#'
#' Table of known repeat expansion disorder loci.
#'
#' @format An exstra_score object:
#' \describe{
#'   \item{data}{A data.table repeat scores, one row for each read}
#'   \item{samples}{A data.table of input samples}
#'   \item{db}{A data.table of the known loci}
#'   \item{input_type}{Indication that the loci are referred to by name (as opposed to genomic position)}
#' }
#' @source 
#'         Rick M. Tankard, Martin B. Delatycki, Paul J. Lockhart, 
#'         Melanie Bahlo. 
#'         Detecting known repeat expansions with standard protocol next generation 
#'         sequencing, towards developing a single screening test for neurological repeat 
#'         expansion disorders. 
#'         \emph{bioRxiv} 157792; 
#'         doi: \url{https://doi.org/10.1101/157792}
"exstra_wgs_pcr_2"
