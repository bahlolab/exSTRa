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
  samples$sex <- factor(NA, c("male", "female"))
  setkey(samples, sample)
  structure(list(data = data.table(data), db = db, samples = samples), class = c("rep_score_data", "strdata"))
}


plot.rep_score_data <- function(rsc, locus = NULL, sample_col = NULL, refline = TRUE, ylab="Fn(x)", verticals = TRUE,
  pch = 19, xlim, ylim = c(0,1), alpha_control = 0.5, alpha_case = NULL, 
  xlinked = "all", xlinked.safe = TRUE, ...) {
  # Plot ECDFs of rep score data
  # sample_col should be a named vector, sample names as the name and color as the value
  # refline: if TRUE, include reference
  # xlinked: For loci on X chromosome, "all" for all samples, "male" and "female" for only that sex
  if(is.null(locus)) {
    strlocis <- strloci(rsc)
  } else {
    strlocis <- locus
  }
  if(!missing(xlim)) {
    xlim_1 <- xlim
  }
  assert('In plot.rep_score_data(), must have xlinked one of "all", "male", "female" or "both"', xlinked %in% c("all", "male", "female", "both"))
  if(xlinked == "both") {
    xlinked_loop <- c("male", "female")
  } else {
    xlinked_loop <- xlinked
  }
  for(locus.name in strlocis) {
    #strrir.trim <- trim.rep_in_read_data(strrir, trimming)
    for(xlinked in xlinked_loop) {
      main.title <- paste(strloci_text_info(rsc$db, locus.name), "score ECDF")
      if(xlinked != "all" && grepl("X", str_score_fil$db$db[locus.name]$Gene.location)) {
        plot_data <- str_filter_sex(rsc, xlinked, xlinked.safe)$data[locus == locus.name]
        main.title <- paste(main.title, paste0(xlinked, 's'))
      } else {
        plot_data <- rsc$data[locus.name]
      }
      if(missing(xlim)) {
        xlim_1 <- c(0, max(plot_data$mlength))
      }
      plot(NA,
        xlim = xlim_1,
        ylim = ylim,
        main = main.title,
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
      for(samp in c(setdiff(unique(plot_data$sample), names(sample_col)), intersect(unique(plot_data$sample), names(sample_col)))) {
        plot(ecdf(plot_data[sample == samp, rep]), add = T, col = replace(sample_col[samp], is.na(sample_col[samp]), black_trans), verticals = verticals,
          pch = pch, ...)
      }
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




rbind.rep_score_data.list <- function(strscore_list, idcol = "data_group") {
  assert("strscore_list must be a list", inherits(strscore_list, "list"))
  if(length(strscore_list) == 0) {
    stop("List is empty")
  }
  assert("Not all elements are rep score data", is.rep_score_data(strscore_list[[1]]))
  if(length(strscore_list) == 1) {
    return(strscore_list[[1]])
  }
  for(i in seq_along(strscore_list)) {
    assert(paste("Element", i, "is not rep_score_data"), is.rep_score_data(strscore_list[[i]]))
    assert("STR database is of mixed types", strscore_list[[1]]$db$input_type == strscore_list[[i]]$db$input_type)
  }
  
  # Could be written much better, all in one go here instead, rather than recursion
  
  data.new <- rbindlist(lapply(strscore_list, function(x) { x$data }), idcol = idcol)
  db.new.db <- rbindlist(lapply(strscore_list, function(x) { x$db$db }))
  setkey(db.new.db, disease.symbol)
  db.new.db <- unique(db.new.db)
  db.new <- strdb(db.new.db, input_type = strscore_list[[1]]$db$input_type)
  new_strscore <- rep_score_data_new(data.new, db.new)
  new_strscore$samples <- rbindlist(lapply(strscore_list, function(x) { x$samples }), idcol = idcol, fill = TRUE)
  return(new_strscore)
}
