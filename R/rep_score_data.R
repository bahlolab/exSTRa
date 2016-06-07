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
  for(locus.name in strlocis) {
    #strrir.trim <- trim.rep_in_read_data(strrir, trimming)
    plot(NA,
      xlim = c(0, 1),
      ylim = c(0, 1),
      main = paste(strloci_text_info(rsc$db, locus.name), "score ECDF"),
      xlab = "Proportion repeated (x)",
      ylab = "Fn(x)",
      cex.main = 1,
      ...)
    grid(col = "grey80")
    plot_data <- rsc$data[locus.name]
    if(is.null(sample_col)) {
      sample_col = ifelse(rsc$samples$group == 'case', "red", "black")
      names(sample_col) <- rsc$samples$sample
    } 
    for(samp in unique(plot_data$sample)) {
      plot(ecdf(plot_data[sample == samp, prop]), add = T, col = replace(sample_col[samp], is.na(sample_col[samp]), "black"))
    }
  }
}

rep_score_data_ks_tests <- function(rsc, locus = NULL, controls = c("control", "all")) {
  # Performs Kolmogorov-Smirnov Tests on samples, comparing to other samples
  #
  # controls allows either just control samples to be used as the population distribution,
  # or all other samples including designated cases. This makes no difference if there is
  # only a single case. 
  if(!is.rep_score_data(rsc)) {
    stop("rsc is not object of class rep_score_data")
  }
  if(is.null(locus)) {
    strlocis <- strloci(rsc)
  } else {
    strlocis <- locus
  } 
  results <- data.table(
    locus = rep(strlocis, length(rsc$samples[group == 'case', sample])), 
    sample = rep(rsc$samples[group == 'case', sample], each = length(strlocis)), 
    p.value = NA_real_ #,
    #test = list(list())
  )
  setkey(results, locus, sample)
  for(loc in strlocis) {
    loc_scores <- rsc$data[locus == loc]
    for(samp in rsc$samples[group == 'case', sample]) {
      KS <- ks.test(loc_scores[group == "control"]$rep, loc_scores[sample == samp]$rep, 
        alternative = "greater", exact = NULL)
      results[list(loc, samp), p.value := KS$p.value]
      #results[list(loc, samp), test := KS]
    }
  }
  
  return(
    data.frame(results)
  )
}
