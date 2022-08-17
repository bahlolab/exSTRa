# Declaring generic methods for S3 classes

#' Return the loci of the selected object
#' 
#' @param x object of an exSTRa class.
#' @param ... additional arguments to loci
#' @return The loci of \code{x}.
#' @examples 
#' loci(exstra_known)
#' loci(exstra_wgs_pcr_2)
#' @export
loci <- function (x, ...) {
  UseMethod("loci", x)
}

#' Give text info for the locus
#' 
#' Usually used in plot titles.
#' @param x object of class exstra_db
#' @param locus Locus to get the information of
#' @export
loci_text_info <- function (x, locus) {
    UseMethod("loci_text_info", x)
}

#' Give the plot names for given sample names
#' 
#' @param x Object
#' @param names Sample names to get the plot names for.
#' @export
plot_names <- function (x, names) {
  UseMethod("plot_names", x)
}

#' Assign the plot names
#' 
#' @param x Object
#' @param value Value to assign.
#' @export
`plot_names<-` <- function (x, value) {
  UseMethod("`plot_names<-`", x)
}

# make data.table copy() also work properly on our class
#' Copy an object with a data.table.
#' 
#' Copy-on-write doesn't apply to data.table. 
#' Instead, an explicit copy() must be used.
#' 
#' @param x Object to copy.
#' @export
copy <- function (x) {
  UseMethod("copy", x)
}

# copy works as normal
#' @export
copy.default <- data.table::copy

