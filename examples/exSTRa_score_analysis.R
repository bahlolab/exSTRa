# An example of exSTRa usage, for known STR expansion disorder loci

## ---- strexpansion_prepare
library(exSTRa)

knitr::opts_chunk$set(fig.width=11, fig.height=11)

# Read score data and file with loci information
str_score <- read_score (
  file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
  database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),  # for greater control, use object from read_exstra_db() instead
  groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), # the group is the first regular expression (regex) to match
  filter.low.counts = TRUE
)

str_score

## ---- Plotting STR score ----
# Plot all loci:
# plot(str_score)

# restrict to only four interesting loci, for simplicity here:
( str_score_four <- str_score[c("HD", "SCA2", "SCA6", "FRDA")] )

# Plot the HD locus only:
plot(str_score["HD"])

# With custom colours:
plot(str_score, "HD", sample_col = c("WGSrpt_10" = "red", "WGSrpt_12" = "blue"))

# For many loci, plot to a file
# Most options not shown here should be passed onto plot.exstra_score() (equivalent to plot() on str_score)
# Without a file.base or directory, should just plot many ECDFs to the R device
# also can take mfrow or mfcol (not both) to place plots in an nr-by-nc array
par(mfrow = c(2, 2))
plot_multi(str_score[c("HD", "SCA6", "FRDA", "SCA1")], dir = "example_images", 
  prefix = "HiSeqXTen_WGS_PCR_2", plot_types = c(1, 2), alpha_case = 0.2)

# For combining data from multiple samples, but analysed in the Perl scripts separately.
data_08 <- str_score[, "WGSrpt_08"]
data_13 <- str_score[, "WGSrpt_13"]
( combined_data <- rbind_score_list(list(data_08, data_13)) )

## ---- Performing tests for expansions ----
# here, the brackets mean the object is shown
(tsum <- tsum_test(str_score_four))

# Plotting tsum only highlights significant samples
plot(tsum)

# You may fix the colours for each sample, as follows: 
plot_cols <- c(RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Dark2"), "orange", "blue")
par(mfrow = c(1, 1))
names(plot_cols) <- str_score_four$samples[, sample]
pie(rep(1, length(plot_cols)), col = plot_cols, labels = names(plot_cols))
plot_cols

par(mfrow = c(2, 2))
# Bonferroni correction is too severe here, so we use Bonferroni correction only on each 
# locus for the number of samples.
plot(tsum, sample_col = plot_cols, correction = "samples")

# Or Bonferroni correction for the number of loci tested:
plot(tsum, sample_col = plot_cols, correction = "loci")

# Give a table of each sample and locus with the p-value, and if it is significant:
(ps <- p_values(tsum, correction = "samples"))

# this may be acted on directly: 
ps[identity(signif)]
# or with the only.signif option:
p_values(tsum, only.signif = TRUE, correction = "samples")
