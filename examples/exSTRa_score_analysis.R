# An example of exSTRa usage, for known STR expansion disorder loci

## ---- strexpansion_prepare
library(exSTRa)

# data.table() # handy if closer inspection of internal tables is required

# Read xlsx database
str_database <- exstra_db_read("data/repeat_disorders.xlsx") # exstra_db object

str_score <- exstra_score_read (
  file = "data/repeat_scores_at_known_loci.txt", # created by Perl script (TODO: exact name)
  database = strdatabase, # alternatively, this may be the direct file path, if no import options are required (TODO)
  groups.regex = c(case = "", control = "^control"), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
  filter.low.counts = TRUE
)

#### ---- STR score ----
# Plot all loci:
# plot(str_score)

# Plot the HD locus only:
plot(str_score["HD"])

# With custom colours
plot(str_score, "HD", sample_col = c("SCA6-1" = "red"))

# For many loci, plot to a file

# Perform tests



