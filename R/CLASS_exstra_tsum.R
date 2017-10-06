# The exstra_tsum class.
# Contains the results and p-values of the test statistic
# Inherits from a exstra_score object, as it contains the raw information leading to the statistic. 
# 
# Importantly, we do not allow modification of the score data at this stage. (To implement)
# Only subsetting may be done by locus. (TODO: maybe not so restrictive, but see what works best)


#' @import data.table
#' @import stringr
#' @import testit
#' 

#' Check if an object is an instance of the exstra_tsum class.
#' 
#' Checks for the \code{class} attribute only. 
#' Does not check for correctness. 
#' 
#' @export
is.exstra_tsum <- function(x) inherits(x, "exstra_tsum")


#' Create a new exstra_tsum object.
exstra_tsum_new_ <- function(strscore, T, pvals) {
  assert("strscore should be from class exstra_score", is.exstra_score(strscore))
  structure(
    list(
      data = strscore$data, 
      db = strscore$db, 
      input_type = strscore$input_type, 
      samples = strscore$samples,
      T = T,
      p.values = pvals
    ), 
    class = c("exstra_tsum", "exstra_score", "exstra_db"))
}

