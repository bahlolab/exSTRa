# Declaring generics

#' @export
loci <- function (x, ...) {
  UseMethod("loci", x)
}

#' @export
loci_text_info <- function (x, ...) {
    UseMethod("loci_text_info", x)
}

#' @export
plot_names <- function (x, ...) {
  UseMethod("plot_names", x)
}

#' @export
`plot_names<-` <- function (x, ...) {
  UseMethod("`plot_names<-`", x)
}

