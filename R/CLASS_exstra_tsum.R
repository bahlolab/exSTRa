# The exstra_tsum class.
# Contains the results and p-values of the test statistic
# Inherits from a exstra_score object, as it contains the raw information leading to the statistic. 
# 
# Importantly, we do not allow modification of the score data at this stage. (To implement)
# Only subsetting may be done by locus. (TODO: maybe not so restrictive, but see what works best)


#' Check if an object is an instance of the exstra_tsum class.
#' 
#' Checks for the class attribute only. 
#' Does not check for correctness. 
#' 
#' @param x Object to be tested
#' 
#' @return Logical
#' 
#' @import data.table
#' @import stringr
#' @import testit
#' @export
is.exstra_tsum <- function(x) inherits(x, "exstra_tsum")


#' Create a new exstra_tsum object.
#' @keywords internal
exstra_tsum_new_ <- function(strscore, tsum, p.values = NULL, 
    qmats = NULL, xecs = NULL, args = NULL, 
  correction = c("bf", "locus", "uncorrected"),
  alpha = 0.05,
  only.signif = FALSE) {
  testit::assert("strscore should be from class exstra_score", is.exstra_score(strscore))
  
  setkey(tsum, locus, sample)
  if(is.null(p.values)) {
    ts_fake <- structure(
      list(stats = tsum, n_tests = tsum[!is.na(tsum), .N])
      , class = c("exstra_tsum", "exstra_score", "exstra_db")
    )
    stats <- p_values(ts_fake,
      correction = correction,
      alpha = alpha,
      only.signif = only.signif,
      modify = TRUE
    )
  } else {
    ps <- p_values(correction = correction,
      alpha = alpha,
      only.signif = only.signif,
      p.matrix = p.values)
    stats <- merge(tsum, ps, all = TRUE)
  }
  
  out_list <- list(
    data = strscore$data, 
    db = strscore$db, 
    input_type = strscore$input_type, 
    samples = strscore$samples,
    stats = stats,
    args = args,
    n_tests = sum (!is.na (stats$tsum))
  )
  # For old tsum_test_1() function:
  if(!is.null(qmats)) {
    out_list <- c(out_list, list(qmats = qmats))
  }
  if(!is.null(xecs)) {
    out_list <- c(out_list, list(xecs = xecs))
  }
  
  structure(
    out_list, 
    class = c("exstra_tsum", "exstra_score", "exstra_db"))
}


#' @export
print.exstra_tsum <- function(x, ...) {
  cat(class(x)[1], " object with ", 
    dim(x$stats)[1], " T sum statistics ($stats),\n  ",
    ifelse(is.null(x$stats$p.value), "without p-values", "with p-values calculated ($stats)"), ",\n",
    "  over ", dim(x$db)[1], ifelse(dim(x$db)[1] == 1, " locus", " loci"), ". ($db)\n",
    sep = "")
  
  if(! is.null(x$stats$p.value)) {
  # exSTRa T := sum of two sample t-tests
  #        
  #                  Raw           Bonferroni (sic) correction 
  # N_p_0.0001:      8             5
  # N_p_0.001:       4             6
  # N_p_0.01:        20            5
  # N_p_0.05:        45            10
  # N_p_remainder:   280           320
  # data: str_score
  # N_samples = 18
  # N_loci = 21
  # N_statistics = 378 - N_NA = 370
  # trim = 0.2
  # 
  cat("\n    T sum statistics summary:\n")
  cat("    exSTRa T := sum of two sample t-tests\n\n")
  
  cat("Alternative hypotheses: subject sample has a larger allele than background samples.\n\n")
  summary_x <- tsum_p_value_summary(x)
  cat("alpha  Bonferroni unadjusted\n")
  for(i in seq_len(summary_x[, .N])) {
    cat(sprintf("%-6g % 10d % 10d", summary_x[i, 1], as.integer(summary_x[i, 2]), as.integer(summary_x[i, 3])), "\n")
  }
  
  # data: str_score
  cat("\n")
  cat("Number of samples:", x$samples[, .N], "\n")
  cat("Number of loci:   ", x$db[, .N], "\n")
  cat("Defined p-values: ", sum(!is.na(x$stats$p.value)), "\n")
  cat("NA p-values:      ", sum(is.na(x$stats$p.value)), "\n")
  if(x$n_tests != sum(!is.na(x$stats$p.value))) {
    cat("Pre-extracted tests:", x$n_tests, "\n")
  }
  cat("Function arguments: trim = ", x$args$trim, 
    ", min.quant = ", x$args$min.quant, 
    ", B = ", x$args$B, 
    "\n", sep = "")
  }
}

#' Plot exstra_tsum, highlighting significant results
#' 
#' Plot the significant results of a tsum
#' 
#' @param x An exstra_tsum object.
#' @param loci Character vector of locus or loci to plot.
#' @param sample_col Named (samples) vector of charcters defining colors. 
#' @param correction Apply this correction method from \code{\link{p_values}}. 
#'            If NULL, then the correction already applied to tsum is used. 
#' @param alpha Significance level. 
#'            If NULL, then the significance already applied to tsum is used. (default 0.05)
#' @param controls_label Generic label for all control samples.
#' @param alpha_nonsignif Transparency alpha value for nonsignificant samples.
#' @param ... further arguments to \code{\link{plot.exstra_score}}
#' 
#' @seealso \code{\link{plot.exstra_score}}
#' 
#' @import data.table
#' @import stringr
#' @import testit
#' @export
plot.exstra_tsum <- function(x, loci = NULL, sample_col = NULL, 
  correction = NULL, alpha = NULL, # when NULL, use significance as-is
  controls_label = "Not significant", 
  alpha_nonsignif = 0.25, 
  ...) {
  # check input
  testit::assert("x should be an exstra_tsum object.", is.exstra_tsum(x))
  if(!is.null(loci)) {
    testit::assert("loci should be a character vector", is.vector(loci), is.character(loci))
    # only work on loci we want to plot
    x <- x[loci]
  }
  if(!is.null(sample_col)) {
    testit::assert("sample_col should be a character vector", is.vector(sample_col), 
      is.character(sample_col))
    testit::assert("sample_col should be named", !is.null(names(sample_col)))
  }
  
  # construct colours
  if(is.null(correction) && is.null(alpha)) {
    # Use the samples marked as significant
    ps <- x$stats[identity(signif)]
  } else {
    if(is.null(correction)) {
      correction <- "bonferroni"
    }
    if(is.null(alpha)) {
      alpha <- 0.05
    }
    ps <- p_values(x, correction = correction, alpha = alpha, only.signif = TRUE)
  }
  significant_sample_colours <- list()
  for(loc in loci(x)) {
    # TODO
    this.ps <- ps[loc] 
    if(is.null(sample_col)) {
      if(this.ps[, .N] > 8) {
        warning("More than 8 significant samples for locus ", loc, 
          ". It may be hard to distiguish samples.")
        this.sample_col <- rainbow(this.ps[, .N])
      } else {
        # we have max() here as brewer.pal() requires at least 3
        this.sample_col <- RColorBrewer::brewer.pal(max(this.ps[, .N], 3), "Set2")
        this.sample_col <- this.sample_col[seq_len(this.ps[, .N])] # for when < 3 signficant samples
      }
      names(this.sample_col) <- this.ps[, sample]
    } else {
      # predefined colours
      this.sample_col <- sample_col[this.ps[, sample]]
    }
    significant_sample_colours[[loc]] <- this.sample_col
  }
  
  # Do the plot:
  # TODO: as.exstra_score(x)
  plot_many_str_score(as.exstra_score(x), loci = loci, 
    plot_cols = significant_sample_colours, 
    controls_label = controls_label, 
    alpha_control = alpha_nonsignif, ...)
  # legend:
  # TODO
}

#' Extract loci or samples from exstra_tsum object
#' 
#' Using \code{i} (select) syntax of data.table to extract loci and/or samples.
#' The first index is the loci filter on x$db, and second sample filter on x$samples. 
#' 
#' @param x exstra_score object
#' @param loc Select loci, using data.table filtering on x$db.
#' @param samp Select samples, using data.table filtering on x$samples.
#' 
#' @return exstra_tsum object
#' 
#' @examples 
#' # Run tsum_test()
#' (tsum <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA2", "SCA6", "FRDA")], B = 100))
#' 
#' # All data:
#' tsum
#' 
#' # Extract one locus:
#' tsum["HD"]
#' 
#' # Extract one sample:
#' tsum[, "WGSrpt_10"]
#' 
#' # One sample at one locus:
#' tsum["HD", "WGSrpt_10"]
#' 
#' # Extract by index 
#' tsum[2, 5]
#' 
#' # The following are intended to work, but do not at present:
#' # Extract all triplet repeats
#' ## tsum[unit_length == 3]
#' ## tsum[unit_length == 3]$db$locus
#' 
#' # Extract all coding repeats
#' ## tsum[gene_region == "coding"]
#' ## tsum[gene_region == "coding"]$db$locus
#' 
#' # Extract all case samples:
#' ## tsum[, group == "case"]
#' 
#' @export
`[.exstra_tsum` <- function(x, loc, samp) {
  testit::assert("locus is not the key of x$db", key(x$db)[1] == "locus")
  # recycle code for exstra_score:
  if(missing(loc)) {
    if(!missing(samp)) {
      x <- `[.exstra_score`(x, , eval(substitute(samp)))
    }
  } else {
    if(missing(samp)) {
      x <- `[.exstra_score`(x, eval(substitute(loc)))
    } else {
      x <- `[.exstra_score`(x, eval(substitute(loc)), eval(substitute(samp)))
    }
  }
  # cut class specific loci
  x$stats <- x$stats[x$db$locus, nomatch=0][sample %in% x$samples$sample, nomatch=0]
  # matrix cut. We do not attempt to filter samples here as it is more complicated, and
  # this is for diagnostics mostly. 
  x$qmats <- x$qmats[x$db$locus]
  x$xecs <- x$xecs[x$db$locus]
  x
}


# Copy an exstra_tsum object
# 
# Prevents changing both objects on changes by reference, that do not copy on write. 
# @param x exstra_tsum object to copy.
# 
# @export
copy.exstra_tsum <- function(x) {
  x <- copy.exstra_score(x)
  x$stats %<>% copy()
  # The following lines probably do nothing, but here in case we change implementation 
  # later on. 
  x$qmats %<>% copy()
  x$xecs %<>% copy()
  x$args %<>% copy()
  x$p.value %<>% copy()
  x
}
