# The rep_score_data class
# Holds scores that give the proportion of a read that matches a given repeat

library(data.table)
library(testit)
library(reshape2)

is.rep_score_data <- function(x) inherits(x, "rep_score_data")

rep_score_data_read <- function(file, database, groups.regex = NULL, groups.samples = NULL) {
  # Load the STR counts
  # Groups should be named null, control and case
  out <- strs_read_(file, database, groups.regex, groups.samples, this.class = "strdata")
  out$data$prop <- with(out$data, rep / mlength)
  return(rep_score_data_new(out$data, out$db))
}


rep_score_data_new <- function(data, db) {
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
  structure(list(data = data.table(data), db = db, samples = samples), class = c("rep_score_data", "strdata"))
}


plot.rep_score_data <- function(rsc, locus = NULL, sample_col = NULL, ...) {
  # Plot ECDFs of rep score data
  # sample_col should be a named vector, sample names as the name and color as the value
  if(is.null(locus)) {
    strlocis <- strloci(rsc)
  } else {
    strlocis <- locus
  }
  show(strlocis)
  for(locus.name in strlocis) {
    #strrir.trim <- trim.rep_in_read_data(strrir, trimming)
    plot(NA,
      xlim = c(0, 1),
      ylim = c(0, 1),
      main = paste(locus.name, "score ECDF"),
      #xlab = "",
      #ylab = "",
      ...)
    grid(col = "grey80")
    plot_data <- rsc$data[locus.name]
    if(is.null(sample_col)) {
      sample_col = ifelse(rsc$samples$group == 'case', "red", "black")
      names(sample_col) <- rsc$samples$sample
    } 
    for(samp in unique(plot_data$sample)) {
      plot(ecdf(plot_data[sample == samp, prop]), add = T, col = sample_col[samp])
    }
  }
}
