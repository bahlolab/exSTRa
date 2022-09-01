#' Give read scores of BAM files
#' 
#' Requires the RSamtools package from Bioconductor to be installed. 
#' (May not work on M1 Macs)
#' 
#' @param paths Paths to BAM files
#' @param sample_names Sample names in the same order as paths. May be set to NULL, 
#'        in which case sample names are found bound sample_name_origin.
#' @param sample_name_origin Option to determine sample names from, either the BAM read group (RG) tag, or path basename. Ignored if sample_names is not NULL.
#' @param sample_name_remove If sample_names == NULL, regular expression passed to stringr::str_remove to 
#'         remove parts of sample names. 
#'         ".bam" extensions are automatically removed. 
#' @param sample_name_extract If sample_names == NULL, regular expression passed to stringr::str_extract to only use part of sample names. 
#' @param scan_bam_flag Sets read filters based on SAM flags. If not NULL, an object returned by Rsamtools::scanBamFlag().
#' @param qname If TRUE, the query name of reads is given in output. 
#' @param method Read score method. 
#'        "overlap" tags bases that form part of a repeat motif (starting on any base), then counts the bases that are tagged with any starting base of the motif. 
#'        "count" simply counts the number of matches of the motif, cycled on its starting bases. 
#'        The original method of the original exSTRa publication closely resembles "overlap", 
#'        though some differences from Bio::STR::exSTRa is expected to occur. 
#'        The "count" method is simplier and therefore probably faster, and may possibly be better. 
#'        We are still evaluating which method is superior. 
#' @param verbosity Control amount of messages, an interger up to 2. 
#' @inheritParams read_score
#' @export
score_bam <- function(paths, database, sample_names = NULL, 
                      sample_name_origin = c("RG", "basename"),
                      sample_name_remove = "",
                      sample_name_extract = ".+",
                      groups.regex = NULL, groups.samples = NULL, 
                      filter.low.counts = TRUE,
                      scan_bam_flag = scanBamFlag(
                        isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
                        isNotPassingQualityControls = FALSE, isDuplicate = FALSE
                      ), 
                      qname = FALSE, method = c("overlap", "count"), 
                      verbosity = 1
              ) {
  if (!require("Rsamtools", quietly = TRUE))
    stop("The package 'Rsamtools' from Bioconductor is required to run this function.")
  method <- match.arg(method)
  sample_name_origin <- match.arg(sample_name_origin)
  checkmate::assert_flag(qname)
  checkmate::assert_flag(filter.low.counts)
  if(!is.null(groups.regex)) checkmate::assert_character(groups.regex)
  #assert("Need groups.samples or groups.regex to be defined", !is.null(groups.samples) || !is.null(groups.regex))
  #assert("Require exactly one of groups.samples or groups.regex to be defined", xor(is.null(groups.samples), is.null(groups.regex)))
  
  if(is.character(database)) {
    checkmate::assert_character(database, len = 1)
    # as database is presumably a file, try to read from it
    database <- read_exstra_db(database)
  } else if(!is.exstra_db(database)) {
    stop("Database is not of a recognised type.")
  }
  
  if(!is.null(sample_names) && length(paths) != length(sample_names)) {
    stop("Length of 'sample_names' does not match length of 'paths'.")
  }
  
  out_list <- list()
  i <- 1
  for(bam_file in paths) {
    if(is.null(sample_names)) {
      if(sample_name_origin == "RG") {
        file_header <- scanBamHeader(bam_file)
        rg <- file_header[[1]]$text[["@RG"]]
        sn <- str_remove(rg[str_detect(rg, "^SM:")], "^SM:")
      } else if(sample_name_origin == "basename") {
        sn <- str_remove(basename(bam_file), "\\.bam$")
      }
      sn <- str_remove(sn, sample_name_remove)
      sn <- str_extract(sn, sample_name_extract)
    } else {
      sn <- sample_names[i]
    }
    if(verbosity >= 1) message("Reading sample ", sn)
    out_list[[sn]] <- score_bam_1(bam_file, database, 
                                 scan_bam_flag = scan_bam_flag, qname = qname,
                                 method = method, verbosity = verbosity)
    i <- i + 1
  }
  counts <- rbindlist(out_list, idcol = "sample")
  # Make like a standard strscore object
  setnames(counts, "qwidth", "mlength")
  
  counts$group <- strs_read_groups_(counts, groups.regex, groups.samples)
  
  strscore <- exstra_score_new_(counts, database)
  
  if(filter.low.counts) {
    # Filter low counts, assumed wanted by default
    strscore <- filter_low_scores(strscore)
  }
  # Remove any loci without data:
  strscore$db <- strscore$db[locus %in% strscore$data$locus]
  
  strscore
}

# Score a single BAM file
score_bam_1 <- function(path, database, sample_names = NULL,
                        scan_bam_flag, qname = FALSE,
                        method = "overlap",
                        verbosity = 1) {
  if(qname) {
    qname_what <- "qname"
  } else {
    qname_what <- character()
  }
  what <- c("rname", "strand", "pos", "qwidth", "seq", "flag", "mapq", qname_what)
  
  bam_head <- scanBamHeader(path)
  
  list_bam_dt <- list()
  for(loc in database$db$locus) {
    if(verbosity >= 2) message("  On locus ", loc)
    whichlist <- list()
    whichlist[[database$db[loc, chrom]]] <- IRanges(database$db[loc, chromStart], database$db[loc, chromEnd])
    which <- do.call(IRangesList, whichlist)
    param <- ScanBamParam(flag = scan_bam_flag, which = which, what = what)
    bam <- scanBam(path, param = param)
    
    bamlist <- list(bam)
    
    # Store names of BAM fields
    bam_field <- names(bamlist[[1]])
    # Go through each BAM field and unlist
    list_loci <- lapply(bam_field, function(y) .unlist(lapply(bamlist, "[[", y))) 
    # Store as data.table 
    bam_dt_list <- purrr::map(list_loci, ~ as.data.table(do.call("DataFrame", .x)))
    bam_dt <- bam_dt_list[[1]]
    # names(bam_dt) <- database$db$locus
    
    motif <- database$db[loc, motif]
    if(database$db[loc, strand == "-"]) {
      motif <- reverseComplement(DNAString(motif)) %>% as.character() # Way to do without conversion that's faster?
    }
    if(method == "overlap") {
      seqmask <- matrix("", nrow = length(bam_dt$seq), ncol = nchar(motif))
      mask <- str_replace_all(motif, ".", ".")
      mcv <- motif_cycles(motif)
      for(jj in seq_along(mcv)) {
        seqmask[, jj] <- str_replace_all(bam_dt$seq, mcv[jj], mask)
      }
      lscore <- rep(0, nrow(seqmask))
      for(ii in seq_len(nrow(seqmask))) {
        # merge the masks
        mask_mat <- str_split(seqmask[ii, ], "", simplify = TRUE)
        mask_seq <- apply(mask_mat, 2, function(x) { ifelse("." %in% x, ".", x[1]) })
        lscore[ii] <- sum(mask_seq == ".")
      }
      bam_dt[, rep := lscore]
      bam_dt[, prop := rep / qwidth]
    } else if(method == "count") {
      lscore <- 0
      for(mc in motif_cycles(motif)) {
        lscore <- lscore + str_count(bam_dt$seq, mc)
      }
      bam_dt[, rep := lscore]
      bam_dt[, prop := rep / (qwidth - nchar(motif) + 1)]
    } else {
      stop("Unknown method. (bug)")
    }
    
    list_bam_dt[[loc]] <- bam_dt
  }
   
  output_table <- rbindlist(list_bam_dt, idcol = "locus")
  output_table[, seq := NULL]
  
  output_table
}

# Function for collapsing the list of lists into a single list
# As per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

# Give all the starting cycles of the motif
motif_cycles <- function(motif) {
  cycles <- character(nchar(motif))
  for(i in seq_along(cycles)) {
    cycles[i] <- paste0(substring(motif, i), substring(motif, 1, i - 1))
  }
  cycles
}
