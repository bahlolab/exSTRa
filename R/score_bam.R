#' Give read scores of BAM files
#' 
#' Requires the RSamtools package from Bioconductor to be installed. 
#' (May not work on M1 Macs)
#' 
#' @param paths Paths to BAM files
#' @param scan_bam_flag Sets read filters based on SAM flags. If not NULL, an object returned by Rsamtools::scanBamFlag().
#' @param qname If TRUE, the query name of reads is given in output. 
#' @inheritParams read_score
#' @export
score_bam <- function(paths, database, sample_names = NULL, 
                      groups.regex = NULL, groups.samples = NULL, 
                      filter.low.counts = TRUE,
                      scan_bam_flag = NULL, qname = FALSE) {
  if (!require("Rsamtools", quietly = TRUE))
    stop("The package 'Rsamtools' from Bioconductor is required to run this function.")
  checkmate::assert_flag(qname)
  checkmate::assert_flag(filter.low.counts)
  if(!is.null(groups.regex)) checkmate::assert_character(groups.regex)
  
  if(is.null(scan_bam_flag)) {
    scan_bam_flag <- scanBamFlag(
      isUnmappedQuery = FALSE, 
      isSecondaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE,
      isDuplicate = FALSE
    )
  }
  out_list <- list()
  out_headers <- list()
  i <- 1
  for(bam_file in paths) {
    out_headers[[i]] <- scanBamHeader(bam_file)
    out_list[[i]] <- score_bam_1(bam_file, database, 
                                 filter.low.counts = filter.low.counts, 
                                 scan_bam_flag = scan_bam_flag, qname = qname)
    i <- i + 1
  }
  list(headers = out_headers, scores = out_list)
}

# Score a single BAM file
score_bam_1 <- function(path, database, sample_names = NULL,
                        filter.low.counts = TRUE,
                        scan_bam_flag, qname = FALSE) {
  
  which <- IRangesList(chr22=IRanges(c(894), c(910)))
  
  which <- GRanges(seqnames = database$db$chrom, IRanges(database$db$chromStart, database$db$chromEnd))
  
  if(qname) {
    qname_what <- "qname"
  } else {
    qname_what <- character()
  }
  what <- c("rname", "strand", "pos", "qwidth", "seq", "flag", "mapq", qname_what)
  param <- ScanBamParam(flag = scan_bam_flag, which = which, what = what)
  
  bam_head <- scanBamHeader(path)
  bam <- scanBam(path, param = param)
  
  bamlist <- list(bam)

  #store names of BAM fields
  bam_field <- names(bamlist[[1]])
  #go through each BAM field and unlist
  list_loci <- lapply(bam_field, function(y) .unlist(lapply(bamlist, "[[", y))) 
  #store as data.table 
  bam_dt <- purrr::map(list_loci, ~ as.data.table(do.call("DataFrame", .x)))
  names(bam_dt) <- exstra_known[c("HD", "SCA1")]$db$locus
  # 
  output_table <- rbindlist(bam_dt, idcol = "locus")
  
  # output_table <- NULL
  
  if(filter.low.counts) {
    # filter_low_scores(str_score)
  }
  list(bam = bam, output_table = output_table)
  output_table
}

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}
