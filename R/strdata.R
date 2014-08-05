# The strdata class
# Holds information related to STR counts as obtained from other programs

library(data.table)
library(testit)

is.strdata <- function(x) inherits(x, "strdata")

strs_read <- function(file, database, groups.regex = NULL, groups = NULL) {
  # Load the STR counts
  # Groups should be named null, control and case
  assert("read.strs requires database to be a strdb", inherits(database, "strdb"))
  assert("Need groups or groups.regex to be defined", !is.null(groups) || !is.null(groups.regex))
  assert("Require exactly one of groups or groups.regex to be defined", xor(is.null(groups), is.null(groups.regex)))
  # Read in data
  counts <- read.delim(file)
  # Add some info to the data
  if(!is.null(groups.regex)) {
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
  data <- data.table(data)
  setkey(data, disease, sample)
  samples <- unique(data[, c("sample", "group"), with = F])
  setkey(samples, sample)
  structure(list(data = data.table(data), db = db, samples = samples), class = c("strdata"))
}

print.strdata <- function(x, ...) {
  cat("strdata object with ", dim(x$data)[1], " observations of type ",  x$db$input_type, ".\n",
      "Includes associated STR database of ", dim(x$db$db)[1], " loci.", sep = "")
}

strloci.strdata <- function(data) {  
  strloci(data$db)
}


