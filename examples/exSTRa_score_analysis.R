# An example of exSTRa usage, for known STR expansion disorder loci

## ---- strexpansion_prepare
library(exSTRa)

# data.table() # handy if closer inspection of internal tables is required

# Read score data and file with loci information
str_score <- exstra_score_read (
  file = "data/HiSeqXTen_WGS_PCR_2.txt", # created by Perl script (TODO: exact name)
  #database = "data/repeat_disorders.xlsx", # for more control, use object from exstra_db_read() instead
  database = "../disease_repeats/repeat_disorders_2017_04_26.xlsx",
  groups.regex = c(case = "^WGSrpt", control = "^WGSrpt_0[24]$"), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
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

## ---- Performing tests for expansions ----


