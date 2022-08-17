#' Return a data.table with p-values of a tsum_exstra object.
#' 
#' 
#' 
#' @param tsum An exstra_tsum object. 
#' @param correction Correction method to use. Use "bf" or TRUE for Bonferroni correction, and 
#'                   "uncorrected" or FALSE for no correction. ("bonferroni" is also acceptable).
#'                   "samples" is Bonferroni correction by the number of tests (samples) at each locus.
#'                   "loci" is Bonferroni correction by the number of loci.
#'                   
#' @param alpha Significance level alpha.
#' @param only.signif If TRUE, only return significant results.
#' @param modify If TRUE, will modify the tsum$stats table. Effectively ignored if only.sig == TRUE.
#' @param p.matrix Matrix of p-values for internal use. Should only be used without tsum. 
#' @return A \code{data.table} keyed by "locus" then "sample". 
#'         Other columns are \code{tsum} as calculated by \code{\link{tsum_test}}, \code{p.value} (uncorrected),
#'         \code{signif} (TRUE if significant after given correction), and 
#'         \code{p.value.sd}, giving the standard deviation of the p-value estimate from the 
#'         simulation.
#' 
#' @import magrittr
#' @import testit
#' @export
p_values <- function(
  tsum, 
  correction = c("bf", "loci", "samples", "uncorrected"),
  alpha = 0.05,
  only.signif = FALSE,
  modify = FALSE, 
  p.matrix = NULL
  ) {
  # verify input
  if(!missing(tsum)) {
    testit::assert("tsum should be an exstra_tsum object.", is.exstra_tsum(tsum))
  }
  testit::assert("correction should be a character or logical vector", 
    is.character(correction) || is.logical(correction), 
    is.vector(correction))
  testit::assert("alpha should be a probability value.", alpha >= 0, alpha <= 1)
  testit::assert("only.signif should be a logical.", is.logical(only.signif))
  testit::assert("", xor(is.null(p.matrix), missing(tsum)))
  
  # Get p.matrix
  if(is.null(p.matrix)) {
    n_tests <- tsum$n_tests
    out.table <- tsum$stats
    if(!modify) { # make a copy so we do not modify it
      out.table %<>% copy()
    }
  } else {
    testit::assert("p.matrix is not a matrix", is.matrix(p.matrix))
    n_tests <- sum (!is.na(p.matrix))
    out.table <- p.matrix %>% 
      reshape2:::melt.matrix(value.name = "p.value", 
                             varnames = c("sample", "locus"), 
                             as.is = TRUE
        ) %>%
      data.table()
  }
  
  # output
  if((is.logical(correction[1]) && correction[1]) || correction[1] == "bf" || correction[1] == "bonferroni") {
    # Bonferroni correction.
    # only correct for tests we have performed:
    out.table[, signif := p.value <= alpha / n_tests ]
  } else if ((is.logical(correction[1]) && ! correction[1]) || correction[1] == "uncorrected") {
    out.table[, signif := p.value <= alpha ]
  } else if (correction[1] == "samples") {
    out.table[!is.na(p.value), N := .N, by = locus]
    out.table[, signif := p.value <= alpha / N ]
    out.table[, N := NULL]
  } else if (correction[1] == "loci") {
    out.table[!is.na(p.value), N := .N, by = sample]
    out.table[, signif := p.value <= alpha / N ]
    out.table[, N := NULL]
  } else {
    stop("Unknown correction method ", correction[1])
  }
  if(only.signif) {
    # Only keep significant results
    out.table <- out.table[identity(signif)]
  }
  setkey(out.table, locus, sample)
  return(out.table[])
}
