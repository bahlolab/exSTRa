# The strdata class
# Holds information related to STR counts as obtained from other programs

library(data.table)
library(testit)
library(reshape2)

is.strdata <- function(x) inherits(x, "strdata")

strs_read_ <- function(file, database, groups.regex = NULL, groups.samples = NULL, this.class = NULL) {
  # Load the STR data, and give it the right class
  assert("read.strs requires database to be a strdb", inherits(database, "strdb"))
  assert("Need groups.samples or groups.regex to be defined", !is.null(groups.samples) || !is.null(groups.regex))
  assert("Require exactly one of groups.samples or groups.regex to be defined", xor(is.null(groups.samples), is.null(groups.regex)))
  assert("This function must have a class to return", !is.null(this.class))
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
  return(list(data = counts, db = database))
}

strs_read <- function(file, database, groups.regex = NULL, groups.samples = NULL) {
  # Load the STR counts
  # Groups should be named null, control and case
  out <- strs_read_(file, database, groups.regex, groups.samples, this.class = "strdata")
  return(strdata_new(out$data, out$db))
}

strdata_new <- function(data, db) {
  assert("data must be of classs data.frame", inherits(data, "data.frame"))
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
  samples$plotname <- NA_character_
  setkey(samples, sample)
  structure(list(data = data.table(data), db = db, samples = samples), class = c("strdata"))
}

print.strdata <- function(x, ...) {
  cat(class(x)[1], " object with ", dim(x$data)[1], " observations of type ",  x$db$input_type, ". ($data)\n",
      "Includes associated STR database of ", dim(x$db$db)[1], " loci. ($db)", sep = "")
}

strloci.strdata <- function(data) {  
  strloci(data$db)
}

set_plotnames <- function(data, labels) {
  assert("data must be of class strdata", inherits(data, "strdata"))
  data$samples[names(labels), plotname := labels]
}

boxplot.strdata <- function(strdata, locus, ..., 
  coverage = NULL,
  read.length = NULL,
  cases.known = FALSE,
  up.cols = c("up_00", "up_01", "up_11", "up_02", "up_12", "up_22"),
  plot.cols = c("up_01", "up_11", "up_02", "up_12"),
  case.x.offset = 0.32
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
  
  disease.info <- tryCatch(strdata$db$db[disease.symbol == locus.in],
    error = function (e) { strdata$db$db[locus == locus.in] } )
  rs.len <- with(disease.info, nchar(as.character(Repeat.sequence)))
  
  ylimits <- NULL
  if(!is.null(coverage)) {
    A1 <- with(disease.info, floor(copyNum * rs.len))
    A2 <- with(disease.info, floor(rn.unst.low * rs.len))
    R <- read.length # for convienience
    exp_case <- numeric()
    if(cases.known) {
      exp_case["up_11"] <- coverage / 2 / R * ( (R <= A1) * (A1 - R + 1) + (R <= A2) * (A2 - R + 1) )
      exp_case["up_02"] <- coverage / 2 / R * ( (R >= A1 + 2) * (R - A1 - 1) + (R >= A2 + 2) * (R - A2 - 1) )
      exp_case["up_01"] <- coverage / 2 / R * ( ((R > A1 + 1) * A1 + (R > A2 + 1) * A2) + ((R <= A1 + 1) + (R <= A2 + 1) * (R - 1)) ) 
    }
    exp_control <- numeric()
    exp_control["up_11"] <- coverage / R * ( (R <= A1) * (A1 - R + 1) )
    exp_control["up_02"] <- coverage / R * ( (R >= A1 + 2) * (R - A1 - 1) )
    exp_control["up_01"] <- coverage / R * ( ( (R > A1 + 1) * A1 ) + (R <= A1 + 1) * (R - 1) )
    ylimits = c(0, max(exp_case, exp_control, features$count))
  }
  
  boxplot(count ~ bin, features[group == "control"], boxwex = 0.4, xaxt='n', 
    xlim = c(.7, 4.4), ylim = ylimits) # , border = "blue")
  axis(1, at = 1:4 + 0.15, labels = levels(features$bin), xlab = "Read location")
  with(features[group == "case"], points(as.numeric(bin) + case.x.offset, count, col = "red"))
  samplenames <- features[group == "case"]$sample
  with(features[group == "case"], text(as.numeric(bin) + case.x.offset - 0.03, count, pos = 4, col = "red", 
    labels = plotnames(strdata, sample)))
  title(strloci_text_info(strdata, locus.in),
    xlab = "Read Location")
  if(!is.null(coverage)) {
    if(cases.known) {
      lines(c(1.3, 1.5), rep(exp_case["up_01"], 2), col = "red")
      lines(c(2.3, 2.5), rep(exp_case["up_11"], 2), col = "red")
      lines(c(3.3, 3.5), rep(exp_case["up_02"], 2), col = "red")
      lines(c(4.3, 4.5), rep(exp_case["up_01"], 2), col = "red")
    }
    lines(c(0.75, 1.25), rep(exp_control["up_01"], 2), col = "blue")
    lines(c(1.75, 2.25), rep(exp_control["up_11"], 2), col = "blue")
    lines(c(2.75, 3.25), rep(exp_control["up_02"], 2), col = "blue")
    lines(c(3.75, 4.25), rep(exp_control["up_01"], 2), col = "blue")
  }
}


plotnames <- function(strdata, names) {
  # gives the plot names for given sample names
  strdata$samples[as.character(names), plotname]
}