## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----ECDF plot HD--------------------------------------------------------

library(exSTRa)
# Read score data and file with loci information
str_score <- read_score (
  file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
  database = system.file("extdata", "repeat_expansion_disorders.txt", package = "exSTRa"),
  groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), 
  filter.low.counts = TRUE
)


# Plot HD locus, With custom colours:
plot(str_score, "HD", sample_col = c("WGSrpt_10" = "red", "WGSrpt_12" = "blue"))
ex.cs1 <- expression(paste("WGSrpt_10"),paste("WGSrpt_12"))
legend(80, 0.3, ex.cs1, lty = 1:1, col = c("blue","red"),  adj = c(0, 0.6),cex = .6)

## ----ECDF plot SCA1------------------------------------------------------

# Plot SCA1 locus, With custom colours:
plot(str_score, "SCA1", sample_col = c("WGSrpt_14" = "red", "WGSrpt_16" = "blue","WGSrpt_08" = "pink"))
ex.cs1 <- expression(paste("WGSrpt_14"),paste("WGSrpt_16"),paste("WGSrpt_08"))
legend(80, 0.3, ex.cs1, lty = 1:1, col = c("blue","red","pink"),  adj = c(0, 0.6),cex = .6)

