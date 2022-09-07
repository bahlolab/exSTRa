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
#' @inheritParams tsum_test
#' @export
score_bam <- function(paths, 
                      database, 
                      sample_names = NULL, 
                      sample_name_origin = c("RG", "basename"),
                      sample_name_remove = "",
                      sample_name_extract = ".+",
                      groups.regex = NULL, 
                      groups.samples = NULL, 
                      filter.low.counts = TRUE,
                      scan_bam_flag = Rsamtools::scanBamFlag(
                        isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
                        isNotPassingQualityControls = FALSE, isDuplicate = FALSE
                      ), 
                      qname = FALSE, 
                      method = c("overlap", "count"), 
                      parallel = FALSE, # TRUE for cluster
                      cluster_n = NULL, # Cluster size if cluster == NULL. When NULL, #threads / 2 (but always at least 1)
                      cluster = NULL, # As created by the parallel package. If cluster == NULL and parallel == TRUE, then a
                      # PSOCK cluster is automatically created with the parallel package.
                      verbosity = 1
              ) {
  if (!requireNamespace("Rsamtools", quietly = TRUE))
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
  
  # Set or check the number of cores in parallel, when no cluster is specified
  if(inherits(cluster, "cluster")) {
    # as a cluster has been given, assume we actually do want to use the parallel package
    parallel <- TRUE
  } else {
    if(parallel) {
      n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
      if(is.null(cluster_n)) {
        # Set the number of cores, max threads / 2 (but at least 1)
        cluster_n <- max(1, n_cores / 2)
      } else {
        if(cluster_n > n_cores) {
          warn.message <- paste0("More threads have been requested by cluster_n (", cluster_n, 
                                 ") than appears to be available (", 
                                 n_cores, ").")
          message(warn.message)
        }
      }
    }
  }
  if(parallel) {
    # Create a new PSOCKcluster if required
    if(is.null(cluster)) {
      # create the cluster, just once
      cluster <- snow::makeCluster(cluster_n)
      on.exit(snow::stopCluster(cluster), add = TRUE)
    }
    
    # Load required functions onto cluster nodes
    snow::clusterEvalQ(cluster, { 
      requireNamespace("magrittr");
      requireNamespace("exSTRa"); 
      requireNamespace("Rsamtools"); 
      requireNamespace("IRanges");
      requireNamespace("stringr");
    })
    snow::clusterExport(cluster, 
      c(".unlist", "motif_cycles", "score_overlap_method", "score_count_method",
        "score_bam_1_locus"))
  }
  
  sample_names <- set_sample_names_score_bam(sample_names, paths, 
                    sample_name_origin, sample_name_remove, sample_name_extract)
  names(paths) <- sample_names
  
  # Run the scoring
  if(parallel) {
    names(paths) <- sample_names
    out_list <- snow::parLapply(cluster, paths, score_bam_1, database, 
                          scan_bam_flag = scan_bam_flag, qname = qname,
                          method = method, verbosity = verbosity)
  } else {
    out_list <- list()
    for(sn in sample_names) {
      if(verbosity >= 1) message("Reading sample ", sn)
      out_list[[sn]] <- score_bam_1(paths[sn], database, 
                                    scan_bam_flag = scan_bam_flag, qname = qname,
                                    method = method, verbosity = verbosity)
    }
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
  
  list_bam_dt <- list()
  for(loc in database$db$locus) {
    if(verbosity >= 2) message("  On locus ", loc)
    list_bam_dt[[loc]] <- score_bam_1_locus(database, loc, path, scan_bam_flag, which, what, method) 
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


# score_overlap_method(bam_dt$seq, motif)
score_overlap_method <- function(seqs, motif) { 
  seqmask <- matrix("", nrow = length(seqs), ncol = nchar(motif))
  mask <- str_replace_all(motif, ".", ".")
  mcv <- motif_cycles(motif)
  for(jj in seq_along(mcv)) {
    seqmask[, jj] <- str_replace_all(seqs, mcv[jj], mask)
  }
  lscore <- rep(0, nrow(seqmask))
  for(ii in seq_len(nrow(seqmask))) {
    # merge the masks
    mask_mat <- str_split(seqmask[ii, ], "", simplify = TRUE)
    mask_seq <- apply(mask_mat, 2, function(x) { ifelse("." %in% x, ".", x[1]) })
    lscore[ii] <- sum(mask_seq == ".")
  }
  lscore
}
  
# score_count_method(bam_dt$seq, motif)
score_count_method <- function(seqs, motif) {
  lscore <- 0
  for(mc in motif_cycles(motif)) {
    lscore <- lscore + str_count(seqs, mc)
  }
  lscore
}


# score_bam_1_locus()
score_bam_1_locus <- function(database, loc, path, scan_bam_flag, which, what, method) {
  requireNamespace("Rsamtools")
  requireNamespace("IRanges")
  whichlist <- list()
  whichlist[[database$db[loc, chrom]]] <- IRanges(database$db[loc, chromStart] - 10, database$db[loc, chromEnd] + 10)
  which <- do.call(IRangesList, whichlist)
  param <- Rsamtools::ScanBamParam(flag = scan_bam_flag, which = which, what = what)
  bam <- Rsamtools::scanBam(path, param = param)
  bamlist <- list(bam)
  
  # Store names of BAM fields
  bam_field <- names(bamlist[[1]])
  # Go through each BAM field and unlist
  list_loci <- lapply(bam_field, function(y) .unlist(lapply(bamlist, "[[", y))) 
  # Store as data.table 
  bam_dt_list <- purrr::map(list_loci, ~ as.data.table(do.call("DataFrame", .x)))
  bam_dt <- bam_dt_list[[1]]
  
  motif <- database$db[loc, motif]
  if(database$db[loc, strand == "-"]) {
    motif <- reverseComplement(DNAString(motif)) %>% as.character() # Way to do without conversion that's faster?
  }
  if(method == "overlap") {
    bam_dt[, rep := score_overlap_method(bam_dt$seq, motif)]
    bam_dt[, prop := rep / qwidth]
  } else if(method == "count") {
    bam_dt[, rep := score_count_method(bam_dt$seq, motif)]
    bam_dt[, prop := rep / (qwidth - nchar(motif) + 1)]
  } else {
    stop("Unknown method. (bug)")
  }
  bam_dt
}

set_sample_names_score_bam <- function(sample_names, paths, sample_name_origin,
                                       sample_name_remove, sample_name_extract){
  if(is.null(sample_names)) {
    sample_names <- character(length(paths))
    i <- 1
    for(bam_file in paths) {
      if(sample_name_origin == "RG") {
        file_header <- scanBamHeader(bam_file)
        rg <- file_header[[1]]$text[["@RG"]]
        if(is.null(rg)) {
          stop("RG line not found in bam: ", bam_file, "\nTry using sample_name_origin = \"basename\"")
        } else {
          sn <- str_remove(rg[str_detect(rg, "^SM:")], "^SM:")
        }
      } else if(sample_name_origin == "basename") {
          sn <- str_remove(basename(bam_file), "\\.bam$")
      }
      if(sample_name_remove != "") {
        sn <- str_remove(sn, sample_name_remove)
      }
      sn <- str_extract(sn, sample_name_extract)
      sample_names[i] <- sn
      i <- i + 1
    }
  }
  sample_names
}