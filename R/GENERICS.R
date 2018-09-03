# Declaring generic methods for S3 classes

#' Return the loci of the selected object
#' 
#' @param x object of an exSTRa class.
#' @return The loci of \code{x}.
#' @examples 
#' loci(exstra_known)
#' loci(exstra_wgs_pcr_2)
#' @export
loci <- function (x, ...) {
  UseMethod("loci", x)
}

#' @export
loci_text_info <- function (x, locus) {
    UseMethod("loci_text_info", x)
}

#' @export
plot_names <- function (x, names) {
  UseMethod("plot_names", x)
}

#' @export
`plot_names<-` <- function (x, value) {
  UseMethod("`plot_names<-`", x)
}

# make data.table copy() also work properly on our class
#' @export
copy <- function (x) {
  UseMethod("copy", x)
}

# copy works as normal
#' @export
copy.default <- data.table::copy

