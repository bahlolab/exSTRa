#' create summary of significance on a exstra_tsum object
#' 
#' Creates a data.table from the given p-values
#' 
#' @param tsum An exstra_tsum object. 
#' @param p P-value thresholds
#' @param bonferroni If TRUE, give Bonferroni corrected values.
#' @param raw If TRUE, give raw p-values.
#' @param bonferroni.size If not NULL, override the number of tests to correct for in the 
#'            Bonferroni correction with this value.
#' 
#' @return A data.table
#' 
#' @seealso \code{\link{p_values}}
#' 
#' @export
tsum_p_value_summary <- function(tsum, 
  p = c(0.0001, 0.001, 0.01, 0.05), 
  bonferroni = TRUE, 
  raw = TRUE, 
  # bonferroni.by.locus = TRUE, # more complicated to implement
  bonferroni.size = NULL) {
  # Check inputs, maybe...
  assert("tsum should be an exstra_tsum object", is.exstra_tsum(tsum))
  if(sum(tsum$stats$p.value > 1, na.rm = TRUE)) {
    stop("Some p-values appear to be above 1. This may be an exSTRa package bug.")
  }
  if(sum(tsum$stats$p.value < 0, na.rm = TRUE)) {
    stop("Some p-values appear to be below 0. This may be an exSTRa package bug.")
  }
  # 
  output <- data.table(alpha = c(p, 1, NA))
  ps <- c(-0.1, p, 1)
  if(bonferroni) {
    ps.bf <- c(-0.1, p / tsum$n_tests, 1)
    output$bf <- 0L
    tab <- table(.bincode(tsum$stats$p.value, ps.bf), useNA = "always")
    output[as.integer(names(tab)), bf := as.integer(tab)]
    output[.N, bf := as.integer(tab[length(tab)])]
  }
  # if(bonferroni.by.locus) {
  #   bf_denominator <- sum (!is.na (tsum$p.values))
  #   ps.bf <- c(-0.1, p / bf_denominator, 1)
  #   output$bf.locus <- 0L
  #   tab <- table(.bincode(tsum$p.values, ps.bf), useNA = "always")
  #   output[as.integer(names(tab)), bf.locus := as.integer(tab)]
  #   output[.N, bf.locus := as.integer(tab[length(tab)])]
  # }
  if(raw) {
    output$raw <- 0L
    tab <- table(.bincode(tsum$stats$p.value, ps), useNA = "always")
    output[as.integer(names(tab)), raw := as.integer(tab)]
    output[.N, raw := as.integer(tab[length(tab)])]
  }
  return(output[])
}
