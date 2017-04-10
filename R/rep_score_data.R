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


plot.rep_score_data <- function(rsc, locus = NULL, sample_col = NULL, refline = TRUE, ylab="Fn(x)", verticals = FALSE,
     pch = 19, xlim, ylim = c(0,1), alpha_control = 0.5, alpha_case = NULL, ...) {
  # Plot ECDFs of rep score data
  # sample_col should be a named vector, sample names as the name and color as the value
  # refline: if TRUE, include reference
  if(is.null(locus)) {
    strlocis <- strloci(rsc)
  } else {
    strlocis <- locus
  }
  if(!missing(xlim)) {
      xlim_1 <- xlim
  }
  for(locus.name in strlocis) {
    #strrir.trim <- trim.rep_in_read_data(strrir, trimming)
    plot_data <- rsc$data[locus.name]
    if(missing(xlim)) {
      xlim_1 <- c(0, max(plot_data$mlength))
    }
    plot(NA,
      xlim = xlim_1,
      ylim = ylim,
      main = paste(strloci_text_info(rsc$db, locus.name), "score ECDF"),
      xlab = "Repeated bases (x)",
      ylab = ylab,
      cex.main = 1,
      ...)
    grid(col = "grey80")
    if(refline) {
      abline(v = strloci_normal_exp(rsc$db, locus.name), col = c("blue", "red"), lty = 3:4)
    }
    black_trans <- rgb(0, 0, 0, alpha = alpha_control)
    if(is.null(sample_col)) {
      sample_col = ifelse(rsc$samples$group == 'case', rgb(1, 0, 0, alpha_case), black_trans)
      names(sample_col) <- rsc$samples$sample
    } 
    if(!is.null(alpha_case)) {
      sample_col <- add.alpha(sample_col, alpha_case)
    }
    for(samp in unique(plot_data$sample)) {
      plot(ecdf(plot_data[sample == samp, rep]), add = T, col = replace(sample_col[samp], is.na(sample_col[samp]), black_trans), verticals = verticals,
        pch = pch, ...)
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

rsd_filter_lower_than_expected <- function(strscore) {
  strscore$db$db[, unit_length := nchar(as.character(Repeat.sequence))]
  # set score, want to remove scores that are smaller than expected by chance
  strscore$db$db[, min_score := unit_length / 4 ^ unit_length]
  strscore$data <- strscore$data[prop > strscore$db$db[as.character(locus), min_score]]
  strscore
}

# TODO:
#this function is from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
#so will need a rewrite
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}
