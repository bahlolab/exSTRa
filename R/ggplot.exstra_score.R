#' @import ggplot2
#' @include exstra_wgs_pcr_2.R
NULL

#' @export
ggplot.exstra_score <- function(data = NULL, ...) {
  #TODO: verify keys
   ggplot(data, ...)
}

#' @export
ggexstra <- function(x, mapping = aes(x = rep, colour = sample), ...) {
  ggplot(x, mapping = mapping, ...) + stat_ecdf()
}



# ggplot(exstra_wgs_pcr_2, aes(x = rep, group = sample)) + stat_ecdf()

# ggplot(exstra_wgs_pcr_2, aes(x = rep, colour = sample)) + stat_ecdf(aes(colour = sample)) + 
#   facet_wrap(~locus)

# ggexstra(exstra_wgs_pcr_2)
