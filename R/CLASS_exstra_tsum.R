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
exstra_tsum_new_ <- function(strscore, T, p.values = NULL, qmats = NULL, xecs = NULL) {
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
      xecs = xecs
    ), 
    class = c("exstra_tsum", "exstra_score", "exstra_db"))
}


#' @export
print.exstra_tsum <- function(x, ...) {
  cat(class(x)[1], " object with ", 
    dim(x$T)[1], " T sum statistics ($T),\n  ",
    ifelse(is.null(x$p), "without p-values", "with p-values calculated ($p.values)"), ",\n",
    "  over ", dim(x$db)[1], ifelse(dim(x$db)[1] == 1, "locus", "loci"), ". ($db)\n",
    sep = "")
}


#TODO:
# brackets [, ]
#' @export
`[.exstra_tsum` <- function(x, loc, samp) {
  assert("locus is not the key of x$T", key(x$T)[1] == "locus")
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
  x$T <- x$T[x$db$locus][sample %in% x$samples$sample]
  x
}
