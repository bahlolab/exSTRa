# An example of exSTRa usage, for known STR expansion disorder loci

## ---- strexpansion_prepare
library(exSTRa)

# data.table() # handy if closer inspection of internal tables is required

# Read score data and file with loci information
str_score <- read_score (
  file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
  database = system.file("extdata", "repeat_expansion_disorders.txt", package = "exSTRa"),  # for more control, use object from read_exstra_db() instead
  groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), # here, matches on successive patterns override previous matches
  filter.low.counts = TRUE
)

str_score

## ---- Plotting STR score ----
# Plot all loci:
# plot(str_score)

# restrict to only three interesting loci, for simplicity here:
( str_score_three <- str_score[c("HD", "SCA6", "FRDA")] )

# Plot the HD locus only:
plot(str_score["HD"])

# With custom colours:
plot(str_score, "HD", sample_col = c("WGSrpt_10" = "red", "WGSrpt_12" = "blue"))
# Add legend:
# TODO

# For many loci, plot to a file
# Most options not shown here should be passed onto plot.exstra_score() (equivalent to plot() on str_score)
# Without a file.base or directory, should just plot many ECDFs to the R device
# also can take mfrow or mfcol (not both) to place plots in an nr-by-nc array
par(mfrow = c(2,2))
plot_multi(str_score_three, dir = "PLACEHOLDER", file.base = "HiSeqXTen_WGS_PCR_2")

( rbinding <- rbind_score_list(list(str_score[, 5], str_score[, 10])) )

## ---- Performing tests for expansions ----
# here, the brackets mean the object is shown
# We use the parallel package to use threads up to one less than available.
(tsum <- tsum_test(str_score_three, parallel = TRUE))


# summarise each locus?
# ideas for options:
summary(tsum, fdr = 0.05)  # by false discovery rate
summary(tsum, p = 0.05)    # by raw p-value
summary(tsum, p_bf = 0.05) # by p-value with bonferroni correction


# if you have positive controls, then if these are specified in
# str_score$samples$pos_control as:
# "LOC" or "LOC1,LOC2,...": where LOC is the locus of a diagnosed expansion in that patient  
#                           Use a comma delimited list for more than one locus
# str_score$samples$neg_control
# "-": confirmed as a negative control for all loci, or can refer to a specific list for
#     all samples str_score$negative_control_loci
# "LOC" or "LOC1,LOC2,...": where LOC[#] are negatively tested loci
# NA for unknown samples
# ROC curves, AUC:
tsum_performance(tsum) 
# restricting to common SCAs and the similar Freidrich (sic) ataxia, likely to have been tested in SCA patients
tsum_performance(tsum, neg_loci = c( "DRPLA", "SCA1", "SCA2", "SCA3", "SCA6", "SCA7", "SCA17", "FRDA"))

