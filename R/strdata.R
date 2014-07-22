# The strdata class
# Holds information related to STR counts as obtained from other programs

library(data.table)
library(testit)

is.strdata <- function(x) inherits(x, "strdata")

strs_read <- function(file, database, groups = NULL, groups.regex = c("control", "case")) {
  # Load the STR counts
  # read.strs(data = "read_counts.txt", database = strdatabase)
  assert("read.strs requires database to be a strdb", inherits(database, "strdb"))
  # Read in data
  counts <- read.delim(file)
  return(strdata_new(counts, database))
}

strdata_new <- function(data, db) {
  assert("data", inherits(data, "data.frame"))
  assert("db must be of class strdb", inherits(db, "strdb"))
  structure(list(data = data.table(data), db = db), class = c("strdata"))
}

print.strdata <- function(x, ...) {
  cat("strdata object with ", dim(x$data)[1], " observations of type ",  x$db$input_type, ".\n",
      "Includes associated STR database of ", dim(x$db$db)[1], " loci.", sep = "")
}

strloci.strdata <- function(data) {  
  strloci(data$db)
}


