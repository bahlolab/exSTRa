# The strdata class
# Holds information related to STR counts as obtained from other programs

library(data.table)
library(testit)

is.strdata <- function(x) inherits(x, "strdata")

strs_read <- function(file, database, groups.regex = NULL, groups.samples = NULL) {
  # Load the STR counts
  # Groups should be named null, control and case
  assert("read.strs requires database to be a strdb", inherits(database, "strdb"))
  assert("Need groups.samples or groups.regex to be defined", !is.null(groups.samples) || !is.null(groups.regex))
  assert("Require exactly one of groups.samples or groups.regex to be defined", xor(is.null(groups.samples), is.null(groups.regex)))
  # Read in data
  counts <- read.delim(file)
  # Add some info to the data
  if(!is.null(groups.regex)) {
    # using regex for groups
    if(is.null(names(groups.regex))) {
      names(groups.regex) <- groups.regex
    }
    assert("Require groups to only to have names 'case', 'control' and 'null'.", is.element(names(groups.regex), c("case", "control", "null")))
    groups_all <- factor(rep(NA, dim(counts)[1]), levels = names(groups.regex))
    for(group.name in names(groups.regex)) {
      groups_all[grepl(groups.regex[group.name], counts$sample)] <- group.name
    }
  }
  if(!is.null(groups.samples)) {
    # using regex for groups
    assert("groups.samples must be a list if used, with vectors with names of at least one of 'case', 'control' or 'null'.", 
           is.list(groups.samples), 
           length(groups.samples) > 0,
           !is.null(names(groups.samples)),
           is.element(names(groups.samples), c("case", "control", "null"))
      )
    assert("groups.samples does not currently accept multiple of the same names for groups, please put all sample names in the one vector under that name", 
           length(unique(names(groups.samples))) == length(names(groups.samples)) )
    if(length(groups.samples) == 1 && names(groups.samples) == "case") {
      # only cases described, so make other samples controls
      groups_all <- factor(rep("control", dim(counts)[1]), levels = c("case", "control"))
      for(sample in groups.samples$case) {
        groups_all[sample == counts$sample] <- "case"
      }      
    } else {
      groups_all <- factor(rep(NA, dim(counts)[1]), levels = names(groups.regex))
      for(group.name in names(groups.samples)) {
        for(sample in groups.samples[[group.name]]) {
          groups_all[sample == counts$sample] <- group.name
        }      
      }     
    }
  }

  counts$group <- as.factor(groups_all)
  return(strdata_new(counts, database))
}

strdata_new <- function(data, db) {
  assert("data", inherits(data, "data.frame"))
  assert("db must be of class strdb", inherits(db, "strdb"))
  data <- data.table(data)
  if(!is.element("locus", colnames(data))) {
    if(is.element("disease", colnames(data))) {
      setnames(data, "disease", "locus")
    } else if(is.element("STR", colnames(data))) {
      setnames(data, "STR", "locus")
    } else {
      stop("Can't find the column with locus name")
    }
  }
  setkey(data, locus, sample) 
  samples <- unique(data[, c("sample", "group"), with = F])
  samples$sample <- as.character(samples$sample)
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


