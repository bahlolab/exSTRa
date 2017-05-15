# An example of exSTRa usage, for known STR expansion disorder loci

## ---- strexpansion_prepare
library(exSTRa)

# data.table() # handy if closer inspection of internal tables is required

# Read score data and file with loci information
str_score <- exstra_score_read (
  file = "data/repeat_scores_at_known_loci_TEMP.txt", # created by Perl script (TODO: exact name)
  #database = "data/repeat_disorders.xlsx", # for more control, use object from exstra_db_read() instead
  database = "../disease_repeats/repeat_disorders_2017_04_26.xlsx",
  groups.regex = c(control = "", case = "^WGSrpt_1[02]$"), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
  filter.low.counts = TRUE
)

#### ---- STR score ----
# Plot all loci:
# plot(str_score)

# Plot the HD locus only:
plot(str_score["HD"])

# With custom colours:
plot(str_score, "HD", sample_col = c("WGSrpt_10" = "red", "WGSrpt_12" = "blue"))

# For many loci, plot to a file

# Perform tests



