# Declaring generics

loci <- function (x, ...) {
  UseMethod("loci", x)
}

loci_text_info <- function (x, ...) {
    UseMethod("loci_text_info", x)
}

plot_names <- function (x, ...) {
  UseMethod("plot_names", x)
}

`plot_names<-` <- function (x, ...) {
  UseMethod("`plot_names<-`", x)
}

