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
#' @source \url{http://www.biorxiv.org/content/early/2017/06/30/157792}
"exstra_wgs_pcr_2"
