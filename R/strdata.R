# The strdata class
# Holds information related to STR counts as obtained from other programs

library(data.table)
library(testit)

is.strdata <- function(x) inherits(x, "strdata")

strs_read <- function(file, database, groups = "control", groups.regex = NULL) {
  # Load the STR counts
  # read.strs(data = "read_counts.txt", database = strdatabase)
  assert("read.strs requires database to be a strdb", inherits(database, "strdb"))
  assert("Need groups or groups.regex to be defined", !is.null(groups) || !is.null(groups.regex))
  assert("Only one of groups or groups.regex to be defined", xor(is.null(groups), is.null(groups.regex)))
  # Read in data
  counts <- read.delim(file)
  # Add some info to the data
  if(is.null(groups)) {
    # using regex for groups
    if(is.null(names(groups.regex))) {
      names(groups.regex) <- groups.regex
    }
    groups <- factor(rep(NA, dim(counts)[1]), levels = names(groups.regex))
    for(group.name in names(groups.regex)) {
      groups[grepl(groups.regex[group.name], counts$sample)] <- group.name
    }
  }
  counts$group <- as.factor(groups)
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


