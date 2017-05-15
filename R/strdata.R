# The strdata class
# Holds information related to STR counts as obtained from other programs

library(data.table)
library(testit)
library(reshape2)

is.strdata <- function(x) inherits(x, "strdata")


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
  samples$sex <- factor(NA, c("male", "female"))
  setkey(samples, sample)
  structure(list(data = data.table(data), db = db, samples = samples), class = c("strdata", "exstra_score"))
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


