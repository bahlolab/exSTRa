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
#' @import data.table
#' @import stringr
#' @import testit
#' @export
is.exstra_tsum <- function(x) inherits(x, "exstra_tsum")


#' Create a new exstra_tsum object.
exstra_tsum_new_ <- function(strscore, T, p.values = NULL, 
    qmats = NULL, xecs = NULL, args = NULL) {
  assert("strscore should be from class exstra_score", is.exstra_score(strscore))
  setkey(T, locus, sample)
  structure(
    list(
      data = strscore$data, 
      db = strscore$db, 
      input_type = strscore$input_type, 
      samples = strscore$samples,
      T = T,
      p.values = p.values,
      qmats = qmats, 
      xecs = xecs,
      args = args,
      n_tests = sum (!is.na (p.values))
    ), 
    class = c("exstra_tsum", "exstra_score", "exstra_db"))
}


#' @export
print.exstra_tsum <- function(x, ...) {
  cat(class(x)[1], " object with ", 
    dim(x$T)[1], " T sum statistics ($T),\n  ",
    ifelse(is.null(x$p.values), "without p-values", "with p-values calculated ($p.values)"), ",\n",
    "  over ", dim(x$db)[1], ifelse(dim(x$db)[1] == 1, " locus", " loci"), ". ($db)\n",
    sep = "")
  
  if(! is.null(x$p.values)) {
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
  # not.signif <- matrix(TRUE, nrow = dim(x$p.values)[1], ncol = dim(x$p.values)[2])
  # not.signif.bf <- matrix(TRUE, nrow = dim(x$p.values)[1], ncol = dim(x$p.values)[2])
  # for(alpha in c(0.0001, 0.001, 0.01, 0.05)) {
  #   is.sig.at.this.level <- not.signif & (x$p.values <= alpha)
  #   is.sig.at.this.level.bf <- not.signif.bf & (x$p.values<= alpha / sum(!is.na(x$p.values)))
  #   cat(alpha, "    ", sum(is.sig.at.this.level.bf, na.rm = TRUE), "    ", sum(is.sig.at.this.level, na.rm = TRUE), "\n") 
  # }
  summary_x <- tsum_p_value_summary(x)
  cat("alpha  Bonferroni unadjusted\n")
  for(i in seq_len(summary_x[, .N])) {
    cat(sprintf("%-6g % 10d % 10d", summary_x[i, 1], as.integer(summary_x[i, 2]), as.integer(summary_x[i, 3])), "\n")
  }
  
  # data: str_score
  cat("\n")
  cat("Number of samples:", x$samples[, .N], "\n")
  cat("Number of loci:   ", x$db[, .N], "\n")
  cat("Defined p-values: ", sum(!is.na(x$p.values)), "\n")
  cat("NA p-values:      ", sum(is.na(x$p.values)), "\n")
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
#' @param tsum An exstra_tsum object.
#' @param locus Character vector of locus or loci to plot.
#' @param sample_col Named (samples) vector of charcters defining colors. 
#' @param ... further arguments to plot.exstra_score
#' 
#' @seealso plot.exstra_score
#' 
#' @export
plot.exstra_tsum <- function(tsum, locus = NULL, sample_col = NULL, 
  correction = "bf", alpha = 0.05, 
  controls_label = "Not significant", ...) {
  # check input
  assert("tsum should be an exstra_tsum object.", is.exstra_tsum(tsum))
  if(!is.null(locus)) {
    assert("locus should be a character vector", is.vector(locus), is.character(locus))
    # only work on loci we want to plot
    tsum <- tsum[locus]
  }
  if(!is.null(sample_col)) {
    assert("sample_col should be a character vector", is.vector(sample_col), 
      is.character(sample_col))
  }
  
  # construct colours
  ps <- p_values(tsum, correction = correction, alpha = alpha, only.signif = TRUE)
  significant_sample_colours <- list()
  for(loc in loci(tsum)) {
    # TODO
    ps
    if(is.null(sample_col)) {
      
    } else {
      
    }
    significant_sample_colours[[loc]] <- NULL # TODO
  }
  
  # Do the plot:
  # TODO: as.exstra_score(tsum)
  plot_many_str_score(as.exstra_score(tsum), locus = locus, 
    plot_cols = significant_sample_colours, 
    controls_label = controls_label, ...)
  # legend:
  # TODO
}

#TODO:
# brackets [, ]
#' @export
`[.exstra_tsum` <- function(x, loc, samp) {
  assert("locus is not the key of x$T", key(x$T)[1] == "locus")
  # recycle code for exstra_score:
  stop("The square bracket function of exstra_tsum objects is not yet written.")
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
  x$T <- x$T[x$db$locus][sample %in% x$samples$sample]
  x
}


#' Copy an exstra_tsum object
#' 
#' @export
copy.exstra_tsum <- function(x) {
  x <- copy.exstra_score(x)
  x$T %<>% copy()
  x$qmats %<>% copy()
  x$xecs %<>% copy()
  x$args %<>% copy()
  x
}
