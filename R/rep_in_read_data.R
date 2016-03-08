# The rep_in_read_data class
# Holds information related to STR (repeat) locations within a read 

library(data.table)
library(testit)
library(reshape2)

is.rep_in_read_data <- function(x) inherits(x, "rep_in_read_data")

rep_in_read_data_read  <- function(file, database, groups.regex = NULL, groups.samples = NULL) {
  # Load the STR counts
  # Groups should be named null, control and case
  out <- strs_read_(file, database, groups.regex, groups.samples, this.class = "strdata")
  return(rep_in_read_data_new(out$data, out$db))
}

rep_in_read_data_new <- function(data, db) {
  assert("data must be of class data.frame", inherits(data, "data.frame"))
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
  structure(list(data = data.table(data), db = db, samples = samples), class = c("rep_in_read_data", "strdata"))
}


trim.rep_in_read_data <- function(data, trim = NULL, trim_a = trim, trim_b = trim) {
  # perform effective trimming of all read data
  assert("Requires input to be rep_in_read_data", is.rep_in_read_data(data))
  assert("Trim must be a non-negative integer.", is.numeric(trim_a), trim_a >= 0, round(trim_a) == trim_a)
  assert("Trim must be a non-negative integer.", is.numeric(trim_b), trim_b >= 0, round(trim_b) == trim_b)
  output <- data # don't want to modify input
  # do the initial trim
  output$data$a <- output$data$a - trim_a
  output$data$c <- output$data$c - trim_b
  # trim into repeat if nessersary, and set 0 for that trimming
  output$data$b <- output$data$b + pmin(output$data$a, 0) + pmin(output$data$c, 0)
  output$data$a <- pmax(output$data$a, 0)
  output$data$c <- pmax(output$data$c, 0)
  # filter reads that no longer contain any repeat
  # Note that if this filter is removed then the trimming will have to take special care
  # as after the algorithm we may have b < 0
  output$data <- output$data[b >= 0]
  message("Trimming removed read data for ", (dim(data$data)[1] - dim(output$data)[1]), " reads of ", dim(data$data)[1])
  return(output)
}


plot.rep_in_read_data <- function(rird, read_length, locus = NULL, sample_col = NULL, inc_text = FALSE, inc_points = TRUE, ...) {
  #
  # sample_col should be a named vector, sample names as the name and color as the value
  strlocis <- ifelse(is.null(locus), strloci(rird), locus)
  for(locus.name in strlocis) {
    #strrir.trim <- trim.rep_in_read_data(strrir, trimming)
    plot(c(0, 1, 1/2, 0) * read_length / sqrt(3/2), c(0, 0, 1, 0) * read_length, xlim = c(0, read_length / sqrt(3/2)),
      ylim = c(0, read_length),
      main = paste(locus.name),
      type = "l",
      col = "grey60", 
      ...)
    grid(col = "grey80")
    plot_data <- rird$data[locus.name]
    if(is.null(sample_col)) {
      dot_col <- ifelse(rird$samples[as.character(plot_data$sample)]$group == 'case', "red", "black")
    } else {
      dot_col <- rep("black", dim(plot_data)[1])
      for(sname in names(sample_col)) {
        dot_col[which(plot_data$sample == sname)] <- sample_col[sname]
      }
    }
    if(inc_points) {
      with(plot_data, 
        points(jitter((a + .5 * b) / sqrt(3/2), 1.5), jitter(b, 1.5), 
          #col = ifelse(rird$samples[get('sample', 2)] == 'case', "red", 
          col = dot_col,
          ...
        )
      ) 
    }
    if(inc_text) {
      with(plot_data, 
        text((a + .5 * b) / sqrt(3/2), b, 
          #col = ifelse(str_rep_in_read$samples[sample.name]$group == "case", "red", "black"),
          labels = sample, 
          cex = 0.5, 
          col = dot_col,
          ...
        )
      )
    }
  }
}



# set_plotnames <- function(data, labels) {
#   assert("data must be of class strdata", inherits(data, "strdata"))
#   data$samples[names(labels), plotname := labels]
# }

# boxplot.rep_in_read_data <- function(rep_in_read_data, locus, ..., 
#   coverage = NULL,
#   read.length = NULL,
#   cases.known = FALSE,
#   up.cols = c("up_00", "up_01", "up_11", "up_02", "up_12", "up_22"),
#   plot.cols = c("up_01", "up_11", "up_02", "up_12"),
#   case.x.offset = 0.32
# ) {
#   # A boxplot for the rep_in_read_data class
#   
#   if(xor(is.null(coverage), is.null(read.length))) {
#     stop("Must specify either none or both of 'coverage' and 'read.length' in boxplot.rep_in_read_data().")
#   }
#   locus.in <- locus
#   stc <- rep_in_read_data$data[locus]
#   
#   features <- melt(stc, id.vars = c("sample", "locus", "group"), variable.name = "bin", 
#     value.name = "count", measure.vars = up.cols)
#   
#   features <- features[is.element(bin, plot.cols)]
#   features$bin <- factor(features$bin, levels = plot.cols)
#   
#   disease.info <- tryCatch(rep_in_read_data$db$db[disease.symbol == locus.in],
#     error = function (e) { rep_in_read_data$db$db[locus == locus.in] } )
#   rs.len <- with(disease.info, nchar(as.character(Repeat.sequence)))
#   
#   ylimits <- NULL
#   if(!is.null(coverage)) {
#     A1 <- with(disease.info, floor(copyNum * rs.len))
#     A2 <- with(disease.info, floor(rn.unst.low * rs.len))
#     R <- read.length # for convienience
#     exp_case <- numeric()
#     if(cases.known) {
#       exp_case["up_11"] <- coverage / 2 / R * ( (R <= A1) * (A1 - R + 1) + (R <= A2) * (A2 - R + 1) )
#       exp_case["up_02"] <- coverage / 2 / R * ( (R >= A1 + 2) * (R - A1 - 1) + (R >= A2 + 2) * (R - A2 - 1) )
#       exp_case["up_01"] <- coverage / 2 / R * ( ((R > A1 + 1) * A1 + (R > A2 + 1) * A2) + ((R <= A1 + 1) + (R <= A2 + 1) * (R - 1)) ) 
#     }
#     exp_control <- numeric()
#     exp_control["up_11"] <- coverage / R * ( (R <= A1) * (A1 - R + 1) )
#     exp_control["up_02"] <- coverage / R * ( (R >= A1 + 2) * (R - A1 - 1) )
#     exp_control["up_01"] <- coverage / R * ( ( (R > A1 + 1) * A1 ) + (R <= A1 + 1) * (R - 1) )
#     ylimits = c(0, max(exp_case, exp_control, features$count))
#   }
#   
#   boxplot(count ~ bin, features[group == "control"], boxwex = 0.4, xaxt='n', 
#     xlim = c(.7, 4.4), ylim = ylimits) # , border = "blue")
#   axis(1, at = 1:4 + 0.15, labels = levels(features$bin), xlab = "Read location")
#   with(features[group == "case"], points(as.numeric(bin) + case.x.offset, count, col = "red"))
#   samplenames <- features[group == "case"]$sample
#   with(features[group == "case"], text(as.numeric(bin) + case.x.offset - 0.03, count, pos = 4, col = "red", 
#     labels = plotnames(rep_in_read_data, sample)))
#   title(strloci_text_info(rep_in_read_data, locus.in),
#     xlab = "Read Location")
#   if(!is.null(coverage)) {
#     if(cases.known) {
#       lines(c(1.3, 1.5), rep(exp_case["up_01"], 2), col = "red")
#       lines(c(2.3, 2.5), rep(exp_case["up_11"], 2), col = "red")
#       lines(c(3.3, 3.5), rep(exp_case["up_02"], 2), col = "red")
#       lines(c(4.3, 4.5), rep(exp_case["up_01"], 2), col = "red")
#     }
#     lines(c(0.75, 1.25), rep(exp_control["up_01"], 2), col = "blue")
#     lines(c(1.75, 2.25), rep(exp_control["up_11"], 2), col = "blue")
#     lines(c(2.75, 3.25), rep(exp_control["up_02"], 2), col = "blue")
#     lines(c(3.75, 4.25), rep(exp_control["up_01"], 2), col = "blue")
#   }
# }


#plotnames <- function(strdata, names) {
#  # gives the plot names for given sample names
#  strdata$samples[as.character(names), plotname]
#}