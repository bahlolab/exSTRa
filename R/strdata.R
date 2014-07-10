# The strdata class
# Holds information related to STR counts as obtained from other programs

library(data.table)
library(testit)

is.strdata <- function(x) inherits(x, "strdata")

strs_read <- function(file, database) {
  # Load the STR counts
  # read.strs(data = "read_counts.txt", database = strdatabase)
  assert("read.strs requires database to be a strdb", inherits(database, "strdb"))
  # Read in data
  counts <- read.delim(data)
  return(strdata_new(counts, database))
}

strdata_new <- function(data, db) {
  structure(list(data = data, db = db), class = c("strdata"))
}
