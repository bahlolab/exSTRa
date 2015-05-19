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

boxplot.strdata <- function(strdata, locus, ..., 
  coverage = NULL,
  read.length = NULL,
  cases.known = FALSE, 
  up.cols = c("up_00", "up_01", "up_11", "up_02", "up_12", "up_22"),
  plot.cols = c("up_01", "up_11", "up_02", "up_12")
) {
  # A boxplot for the strdata class
  
  if(xor(is.null(coverage), is.null(read.length))) {
    stop("Must specify either none or both of 'coverage' and 'read.length' in boxplot.strdata().")
  }
  locus.in <- locus
  stc <- strdata$data[locus]
  
  features <- melt(stc, id.vars = c("sample", "locus", "group"), variable.name = "bin", 
    value.name = "count", measure.vars = up.cols)
  
  features <- features[is.element(bin, plot.cols)]
  features$bin <- factor(features$bin, levels = plot.cols)
  
  disease.info <- strdata$db$db[disease.symbol == locus.in]
  rs.len <- with(disease.info, nchar(as.character(Repeat.sequence)))
  
  ylimits <- NULL
  if(!is.null(coverage)) {
    A1 <- with(disease.info, floor(copyNum * rs.len))
    A2 <- ifelse(cases.known, with(disease.info, floor(rn.unst.low * rs.len)), A1)
    R <- read.length # for convienience 
    ex_up_11 <- coverage / 2 / R * ( (R <= A1) * (A1 - R + 1) + (R <= A2) * (A2 - R + 1) )
    ex_up_02 <- coverage / 2 / R * ( (R >= A1 + 2) * (R - A1 - 1) + (R >= A2 + 2) * (R - A2 - 1) )
    ex_up_01 <- coverage / 2 * ( ( (R > A1 + 1) * A1 + (R > A2 + 1) * A2) / R + (R <= A1 + 1) + (R <= A2 + 1))
    ylimits = c(0, max(ex_up_11, ex_up_02, ex_up_01, features$count))
  }
  
  boxplot(count ~ bin, features[group == "control"], boxwex = 0.5, xaxt='n', 
    xlim = c(.7, 4.4), ylim = ylimits)
  axis(1, at = 1:4 + 0.15, labels = levels(features$bin), xlab = "Read location")
  with(features[group == "case"], points(as.numeric(bin) + 0.34, count, col = "red"))
  title(with(disease.info, 
    paste0(locus.in, " (", 
      Location.of.repeat.within.gene, " ", Repeat.sequence, ") norm: ", floor(copyNum), 
      " (", floor(copyNum * rs.len),"bp) , exp: ", rn.unst.low, " (", 
      floor(rn.unst.low * rs.len), "bp)")),
    xlab = "Read Location")
  if(!is.null(coverage)) {
    lines(c(0.7, 1.5), rep(ex_up_01, 2), col = "blue")
    lines(c(1.7, 2.5), rep(ex_up_11, 2), col = "blue")
    lines(c(2.7, 3.5), rep(ex_up_02, 2), col = "blue")
    lines(c(3.7, 4.5), rep(ex_up_01, 2), col = "blue")
    #TODO: adjust graph so we can see higher values
  }
  
}



