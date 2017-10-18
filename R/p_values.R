#' Return a data.table with p-values of a tsum_exstra object.
#' 
#' 
#' 
#' @param tsum An exstra_tsum object
#' @param correction Correction method to use. Use "bf" or TRUE for Bonferroni correction, and 
#'                   "uncorrected" or FALSE for no correction. ("bonferroni" is also acceptable).
#'                   "locus" is Bonferroni correction by locus.
#'                   
#' @param alpha Significance level alpha.
#' @param only.sig If TRUE, only return significant results.
#' @return A data.table
#' 
#' @import magrittr
#' @import testit
#' @export
p_values <- function(
  tsum, 
  correction = c("bf", "locus", "uncorrected"),
  alpha = 0.05,
  only.signif = FALSE
  ) {
  # verify input
  assert("tsum should be an exstra_tsum object.", is.exstra_tsum(tsum))
  assert("correction should be a character or logical vector", 
    is.character(correction) || is.logical(correction), 
    is.vector(correction))
  assert("alpha should be a probability value.", alpha >= 0, alpha <= 1)
  assert("only.signif should be a logical.", is.logical(only.signif))
  
  # output
  out.table <- tsum$p.values %>% 
    reshape2:::melt.matrix(value.name = "p.value", varnames = c("sample", "locus")) %>%
    data.table()
  if((is.logical(correction[1]) && correction[1]) || correction[1] == "bf" || correction[1] == "bonferroni") {
    # Bonferroni correction.
    # only correct for tests we have performed:
    out.table[, signif := p.value <= alpha / tsum$n_tests ]
  } else if ((is.logical(correction[1]) && ! correction[1]) || correction[1] == "uncorrected") {
    out.table[, signif := p.value <= alpha ]
  } else if (correction[1] == "locus") {
    out.table[!is.na(p.value), N := .N, by = locus]
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
