library(data.table)
library(exSTRa)
library(ggplot2)

ggplot.exstra_score <- function(data = NULL, ...) {
  
   ggplot(data$data, 
    ...)
}

ggplot(exstra_wgs_pcr_2, aes(x = rep, group = sample)) + stat_ecdf()

ggplot(exstra_wgs_pcr_2, aes(x = rep, colour = sample)) + stat_ecdf(aes(colour = sample)) + 
  facet_wrap(~locus)
