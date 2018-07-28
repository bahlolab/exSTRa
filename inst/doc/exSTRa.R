## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6,
  dev = "png",
  fig.retina = 2
)

## ------------------------------------------------------------------------
library(data.table) # should be loaded before exSTRa
library(exSTRa)

## ------------------------------------------------------------------------
str_score <- read_score (
  file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
  database = system.file("extdata", "repeat_expansion_disorders.txt", package = "exSTRa"),
  groups.regex = c(control = "^WGSrpt_0[24]$", case = "")
)

# an exstra_score object:
str_score

## ------------------------------------------------------------------------
str_score$db$locus

## ------------------------------------------------------------------------
plot(str_score["HD"])

## ------------------------------------------------------------------------
plot(str_score, "HD", sample_col = c("WGSrpt_10" = "red", "WGSrpt_12" = "blue"))

## ------------------------------------------------------------------------
( str_score_four <- str_score[c("HD", "SCA2", "SCA6", "FRDA")] )

## ---- out.width = '82%', fig.width=12, fig.height=12---------------------
par(mfrow = c(2, 2))
plot_multi(str_score_four, dir = "example_images", 
  prefix = "HiSeqXTen_WGS_PCR_2", plot_types = c(1, 2), alpha_case = 0.2)

## ------------------------------------------------------------------------
( tsum <- tsum_test(str_score_four, parallel = TRUE) )

## ---- out.width = '82%', fig.width=12, fig.height=12---------------------
par(mfrow = c(2, 2))
plot(tsum)

## ---- fig.width=6, fig.height=6------------------------------------------
plot_cols <- c(RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Dark2"), "orange", "blue")
names(plot_cols) <- str_score_four$samples[, sample]
plot_cols
# To demonstrate the colours used:
par(mfrow = c(1, 1))
pie(rep(1, length(plot_cols)), col = plot_cols, labels = names(plot_cols), cex = 0.7)

## ---- out.width = '82%', fig.width=12, fig.height=12---------------------
par(mfrow = c(2, 2))
plot(tsum, sample_col = plot_cols, correction = "samples")

## ---- out.width = '82%', fig.width=12, fig.height=12---------------------
par(mfrow = c(2, 2))
plot(tsum, sample_col = plot_cols, correction = "loci")

## ------------------------------------------------------------------------
(ps <- p_values(tsum, correction = "samples"))

## ------------------------------------------------------------------------
ps[identity(signif)]

## ------------------------------------------------------------------------
p_values(tsum, only.signif = TRUE, correction = "samples")

## ------------------------------------------------------------------------
exstra_wgs_pcr_2["HD"]

## ------------------------------------------------------------------------
exstra_wgs_pcr_2[, "WGSrpt_10"]

## ------------------------------------------------------------------------
exstra_wgs_pcr_2[unit_length == 3]
exstra_wgs_pcr_2[unit_length == 3]$db$locus

## ------------------------------------------------------------------------
exstra_wgs_pcr_2[gene_region == "coding"]
exstra_wgs_pcr_2[gene_region == "coding"]$db$locus

## ------------------------------------------------------------------------
exstra_wgs_pcr_2[, group == "case"]

## ------------------------------------------------------------------------
data_08 <- str_score[, "WGSrpt_08"]
data_13 <- str_score[, "WGSrpt_13"]
( combined_data <- rbind_score_list(list(data_08, data_13)) )

