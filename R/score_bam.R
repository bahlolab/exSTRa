#' Give read scores of BAM files
#' 
#' Requires the RSamtools package from Bioconductor to be installed. 
#' (May not work on M1 Macs)
#' 
#' @param paths Paths to BAM files
#' @inheritParams read_score
#' @export
score_bam <- function(paths, database, sample_names = NULL, 
                      groups.regex = NULL, groups.samples = NULL, 
                      filter.low.counts = TRUE) {
  if (!require("Rsamtools", quietly = TRUE))
    stop("The package 'Rsamtools' from Bioconductor is required to run this function.")
  
  bam_file <- paths[1]
  
  which <- IRangesList(seq2=IRanges(c(894), c(910)))
  what <- c("rname", "strand", "pos", "qwidth", "seq")
  param <- ScanBamParam(which=which, what=what)

  
  bam <- scanBam(bam_file, param=param)
  
}

