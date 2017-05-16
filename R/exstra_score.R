# The exstra_score class
# Holds scores that give the proportion of a read that matches a given repeat

library(data.table)
library(testit)
# library(reshape2) # I forget where I use this...

is.exstra_score <- function(x) inherits(x, "exstra_score")
# make this the main class

strs_read_ <- function(file, database, groups.regex = NULL, groups.samples = NULL, this.class = NULL) {
  # Load the STR data, and give it the right class
  assert("read.strs requires database to be class exstra_db", is.exstra_db(database))
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

exstra_score_read <- function(file, database, groups.regex = NULL, groups.samples = NULL, filter.low.counts = TRUE) {
  # Load the STR counts
  # Groups should be named null, control and case
  if(is.character(database)) {
    # as database is presumbly a file, try to read from it
    database <- exstra_db_read(database)
  }
  out <- strs_read_(file, database, groups.regex, groups.samples, this.class = "exstra_score")
  # TODO: checks for this input
  out$data$prop <- with(out$data, rep / mlength)
  strscore <- exstra_score_new_(out$data, out$db)
  if(filter.low.counts) {
    # Filter low counts, assumed wanted by default
    strscore <- exstra_low_filter(strscore)
  }
  return(strscore)
}


exstra_score_new_ <- function(data, db) {
  assert("data must be of class data.frame", inherits(data, "data.frame"))
  assert("db must be of class exstra_db", inherits(db, "exstra_db"))
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
  structure(list(data = data.table(data), db = db, samples = samples), class = c("exstra_score"))
}

print.exstra_score <- function(x, ...) {
  cat(class(x)[1], " object with ", dim(x$data)[1], " observations of type ",  x$db$input_type, "($data),\n",
    "  for ", dim(x$samples)[1], " samples. ($samples)\n",
    "  Includes associated STR database of ", dim(x$db$db)[1], " loci. ($db)\n", 
    sep = "")
}

loci.exstra_score <- function(data) {  
  loci(data$db)
}

plot_names.exstra_score <- function(strscore, names) {
  # gives the plot names for given sample names
  strscore$samples[as.character(names), plotname]
}

`plot_names.exstra_score<-` <- function(data, labels) {
  assert("data must be of class exstra_db", inherits(data, "exstra_db"))
  data$samples[names(labels), plotname := labels]
}


# This function allows square brackets to be used to select out the locus and sample
# BIG TODO: always list by locus, throughout the code!!!
`[.exstra_score` <- function(x, loc, samp) {
  assert("locus is not the key of x$data", key(x$data)[1] == "locus")
  assert("sample is not the key of x$samples", key(x$samples)[1] == "sample")
  assert("locus not the key of x$db$db", key(x$db$db)[1] == "locus")
  if(!missing(loc)) {
    x$db$db <- x$db$db[eval(substitute(loc))]
  }
  if(!missing(samp)) {
    x$samples <- x$samples[eval(substitute(samp))]
  }
  x$data <- x$data[x$db$db$locus][sample %in% x$samples$sample]
  x
}

# old name: str_filter_sex
# filter rep_score_data by sex
exstra_filter_sex <- function(strscore, sex = "known", safe = TRUE) {
  # filter rep_score_data by sex
  # sex can be:
  #   "all":     no filtering
  #   "male":    only male samples
  #   "female":  only female samples
  #   "missing": only missing samples
  #   "known":   only samples with sex assigned
  # When safe is TRUE, missing sex assignments with cause an error for sex filtering of 
  #    "all", "male" or "female"
  if(sex %in% c("all", "male", "female")) {
    if(safe) {
      # Check that no data is missing
      if(sum(is.na(strscore$samples$sex)) != 0) {
        stop("In str_filter_sex(), some samples have not been assigned a sex.")
      }
    }
    if(sex == "all") {
      return(strscore)
    } else if(sex == "male") {
      return(strscore[, sex == "male"])
    } else if(sex == "female") {
      return(strscore[, sex == "female"])
    }
  } else if (sex == "missing") {
    strscore[, is.na(sex)]
  } else if (sex == "known") {
    strscore[, !is.na(sex)]
  } else {
    stop("Bad sex assignment")
  }
}

plot.exstra_score <- function(rsc, locus = NULL, sample_col = NULL, refline = TRUE, ylab="Fn(x)", verticals = TRUE,
  pch = 19, xlim, ylim = c(0,1), alpha_control = 0.5, alpha_case = NULL, 
  xlinked = "all", xlinked.safe = TRUE, ...) {
  # Plot ECDFs of rep score data
  # sample_col should be a named vector, sample names as the name and color as the value
  # refline: if TRUE, include reference
  # xlinked: For loci on X chromosome, "all" for all samples, "male" and "female" for only that sex
  if(is.null(locus)) {
    strlocis <- loci(rsc)
  } else {
    strlocis <- locus
  }
  if(!missing(xlim)) {
    xlim_1 <- xlim
  }
  assert('In plot.exstra_score(), must have xlinked one of "all", "male", "female" or "both"', xlinked %in% c("all", "male", "female", "both"))
  if(xlinked == "both") {
    xlinked_loop <- c("male", "female")
  } else {
    xlinked_loop <- xlinked
  }
  for(locus.name in strlocis) {
    #strrir.trim <- trim.rep_in_read_data(strrir, trimming)
    for(xlinked in xlinked_loop) {
      main.title <- paste(loci_text_info(rsc$db, locus.name), "score ECDF")
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

# old name: rsd_filter_lower_than_expected
# Filter read scores with lower than expected scores, under
# the assumption each base in the sequence is uniform and independent.
exstra_low_filter  <- function(strscore) {
  strscore$db$db[, unit_length := nchar(as.character(Repeat.sequence))]
  # set score, want to remove scores that are smaller than expected by chance
  strscore$db$db[, min_score := unit_length / 4 ^ unit_length]
  small_db <- strscore$db$db[, list(locus, min_score)]
  setkey(small_db, locus)
  # strscore$data <- strscore$data[prop > strscore$db$db[as.character(locus), min_score]]
  strscore$data <- strscore$data[small_db][prop > min_score][, min_score := NULL]
  setkey(strscore$data, locus, sample)
  #TODO: check bizarre behaviour of data not printing first time here...
  strscore
}


# TODO:
# this function is from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
# so may need a rewrite
# old name: add.alpha
add_alpha_ <- function(col, alpha = 1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  if(length(col) == 0) {
    stop("col is of length 0")
  }
  apply(sapply(col, col2rgb)/255, 2, 
    function(x) 
      rgb(x[1], x[2], x[3], alpha=alpha))  
}


# Combine multiple exstra_score objects, checking for sample name clashes
# old name: rbind.exstra_score.list
rbind_exstra_score_list <- function(strscore_list, idcol = "data_group", allow_sample_clash = FALSE) {
  assert("strscore_list must be a list", inherits(strscore_list, "list"))
  if(length(strscore_list) == 0) {
    stop("List is empty")
  }
  assert("Not all elements are rep score data", is.exstra_score(strscore_list[[1]]))
  if(length(strscore_list) == 1) {
    return(strscore_list[[1]])
  }
  for(i in seq_along(strscore_list)) {
    assert(paste("Element", i, "is not exstra_score"), is.exstra_score(strscore_list[[i]]))
    assert("STR database is of mixed types", strscore_list[[1]]$db$input_type == strscore_list[[i]]$db$input_type)
  }
  
  # Could be written much better, all in one go here instead, rather than recursion
  
  data.new <- rbindlist(lapply(strscore_list, function(x) { x$data }), idcol = idcol)
  db.new.db <- rbindlist(lapply(strscore_list, function(x) { x$db$db }))
  setkey(db.new.db, locus)
  db.new.db <- unique(db.new.db)
  db.new <- exstra_db(db.new.db, input_type = strscore_list[[1]]$db$input_type)
  new_strscore <- exstra_score_new_(data.new, db.new)
  new_strscore$samples <- rbindlist(lapply(strscore_list, function(x) { x$samples }), idcol = idcol, fill = TRUE)
  setkey(new_strscore$samples, sample)
  if(!allow.sample.clash) {
    test <- table (new_strscore$samples$sample)
    if(max(test) > 1) {
      stop("A sample name is duplicated in inputs, for sample names: ", 
        paste(names(which(test > 1)), collapse = ", "), 
        "\nSet allow.sample.clash = TRUE if this is ok. "
      )
    }
  }
  return(new_strscore)
}

# convinient version of rbind_exstra_score_list() without the use of lists
rbind_exstra_score <- function(..., idcol = "data_group", allow_sample_clash = FALSE) {
  rbind_exstra_score_list(list(...), idcol = idcol, allow_sample_clash = allow_sample_clash)
}

# TODO: easy renaming of samples

# Functions likely of no use:
exstra_score_ks_tests <- function(rsc, locus = NULL, controls = c("control", "all")) {
  # Performs Kolmogorov-Smirnov Tests on samples, comparing to other samples
  #
  # controls allows either just control samples to be used as the population distribution,
  # or all other samples including designated cases. This makes no difference if there is
  # only a single case. 
  if(!is.exstra_score(rsc)) {
    stop("rsc is not object of class exstra_score")
  }
  if(is.null(locus)) {
    strlocis <- loci(rsc)
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


loci_text_info.exstra_score <- function(x, ...) {
  loci_text_info(x$db, ...)
}
