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
  if (!require("GenomicAlignments", quietly = TRUE))
    stop("The package 'GenomicAlignments' from Bioconductor is required to run this function.")
  if (!require("Rsamtools", quietly = TRUE))
    stop("The package 'Rsamtools' from Bioconductor is required to run this function.")
  
  out_list <- list()
  out_headers <- list()
  i <- 1
  for(bam_file in paths) {
    out_headers[[i]] <- scanBamHeader(bam_file)
    out_list[[i]] <- score_bam_1(bam_file, database, filter.low.counts = filter.low.counts)
    i <- i + 1
  }
  list(headers = out_headers, scores = out_list)
}

# Score a single BAM file
score_bam_1 <- function(path, database, sample_names = NULL, 
                        groups.regex = NULL, groups.samples = NULL, 
                        filter.low.counts = TRUE) {
  
  which <- IRangesList(seq2=IRanges(c(894), c(910)))
  what <- c("rname", "strand", "pos", "qwidth", "seq")
  param <- ScanBamParam(which=which, what=what)
  
  bampal <- readGAlignmentsList(path, param = param, with.which_label = TRUE)
  
  
  count_list <- list()
  # worry about optimization later
  for(i in seq_along(bampal)) {
    count_list[[i]] <- str_count(bampal[[i]]@elementMetadata$seq, "CAG")
  }
  str_score <- unlist(count_list)
  
  if(filter.low.counts) {
    # filter_low_scores(str_score)
  }
  str_score
}

