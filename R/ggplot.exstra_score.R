#' Default ggplot method for exstra_score
#'
#' @param data exstra_score object
#' @param ... Additional arguments to ggplot()
#' @examples 
#' library(ggplot2)
#' ggplot(exstra_wgs_pcr_2, aes(x = rep, colour = sample)) +
#'   stat_ecdf() +
#'   facet_wrap("locus")
#' 
#' @include exstra_wgs_pcr_2.R
#' @import ggplot2
#' @export
ggplot.exstra_score <- function(data = NULL, ...) {
  verify.exstra_score(data)
  # Merge sample data into main data
  data <- merge(data$data, data$samples, by = "sample")
  ggplot(data, ...)
}

#' Create a ggplot2 ecdf plot from an exstra_score object
#' 
#' @param x exstra_score object
#' @param mapping Default list of aesthetic mappings to use for plot.
#' @param ... Additional arguments to ggplot()
#' @examples 
#' ggexstra_ecdf(exstra_wgs_pcr_2["HD",])
#' 
#' library(ggplot2)
#' ggexstra_ecdf(exstra_wgs_pcr_2[,]) + facet_wrap("locus")
#'
#' @include exstra_wgs_pcr_2.R
#' @import ggplot2
#' @export
ggexstra_ecdf <- function(x, mapping = aes(x = rep, colour = sample), ...) {
  ggplot(x, mapping = mapping, ...) + stat_ecdf()
}



# ggplot(exstra_wgs_pcr_2, aes(x = rep, group = sample)) + stat_ecdf()

# ggplot(exstra_wgs_pcr_2, aes(x = rep, colour = sample)) + stat_ecdf(aes(colour = sample)) + 
#   facet_wrap(~locus)

# ggexstra(exstra_wgs_pcr_2)
