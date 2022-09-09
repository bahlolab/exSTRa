# The exstra_score class
# Holds scores that give the proportion of a read that matches a given repeat.

#' @import data.table
#' @import stringr
#' @import testit
#' 
#' @title Test if object is an exstra_score object
#' 
#' @param x Object to be tested
#' 
#' @return Logical
#' 
#' @export
is.exstra_score <- function(x) inherits(x, "exstra_score")

#' @importFrom utils read.delim
strs_read_ <- function(file, database, groups.regex = NULL, groups.samples = NULL, this.class = NULL) {
  # Load the STR data, and give it the right class
  assert("read.strs requires database to be class exstra_db", is.exstra_db(database))
  assert("This function must have a class to return", !is.null(this.class))
  # Read in data
  counts <- read.delim(file)
  # Add some info to the data
  counts$group <- strs_read_groups_(counts, groups.regex, groups.samples)
  return(list(data = counts, db = database))
}

# Find groups
strs_read_groups_ <- function(counts, groups.regex = NULL, groups.samples = NULL) {
  assert("Need groups.samples or groups.regex to be defined", !is.null(groups.samples) || !is.null(groups.regex))
  assert("Require exactly one of groups.samples or groups.regex to be defined", xor(is.null(groups.samples), is.null(groups.regex)))
  if(!is.null(groups.regex)) {
    # using regex for groups
    groups.regex <- rev(groups.regex) # we want the first argument of groups.regex to take priority, this behaviour replaces the old behaviour
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
    # Specify the group of samples directly via groups.samples
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
      groups_all <- factor(rep(NA, dim(counts)[1]), levels = names(groups.samples))
      for(group.name in names(groups.samples)) {
        for(sample in groups.samples[[group.name]]) {
          groups_all[sample == counts$sample] <- group.name
        }      
      }     
    }
  }
  return(as.factor(groups_all))
}


#' Create a new exstra_score object
#' @keywords internal
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
  samples <- unique(data[, c("sample", "group"), with = FALSE])
  samples$sample <- as.character(samples$sample)
  samples$plotname <- NA_character_
  samples$sex <- factor(NA, c("male", "female"))
  setkey(samples, sample)
  structure(list(data = data.table(data), db = db$db, input_type = db$input_type, samples = samples), class = c("exstra_score", "exstra_db"))
}

#' @export
print.exstra_score <- function(x, ...) {
  cat(class(x)[1], " object with ", dim(x$data)[1], " observations of type ",  x$input_type, " ($data),\n",
    "  for ", dim(x$samples)[1], " samples. ($samples)\n",
    "  Includes associated STR database of ", 
       dim(x$db)[1], ifelse(dim(x$db)[1] == 1, " locus", " loci"), ". ($db)\n", 
    sep = "")
}

#' @export
plot_names.exstra_score <- function(x, names = NULL) {
  # gives the plot names for given sample names
  if(is.null(names)) {
    x$samples[, plotname]
  } else {
    x$samples[as.character(names), plotname]
  }
}

#' @export
`plot_names<-.exstra_score` <- function(x, value) {
  assert("data must be of class exstra_db", inherits(x, "exstra_db"))
  x$samples[names(value), plotname := value]
}


# This function allows square brackets to be used to select out the locus and sample
# BIG TODO: always list by locus, throughout the code!!!
#' Extract loci or samples frome exstra_score object
#' 
#' Using \code{i} (select) syntax of data.table to extract loci and/or samples.
#' The first index is the loci filter on x$db, and second sample filter on x$samples. 
#' 
#' @param x exstra_score object
#' @param loc Select loci, using data.table filtering on x$db.
#' @param samp Select samples, using data.table filtering on x$samples.
#' 
#' @return exstra_score object
#' 
#' @examples 
#' 
#' # All data:
#' exstra_wgs_pcr_2
#' 
#' # Extract one locus:
#' exstra_wgs_pcr_2["HD"]
#' 
#' # Extract one sample:
#' exstra_wgs_pcr_2[, "WGSrpt_10"]
#' 
#' # Extract all triplet repeats
#' exstra_wgs_pcr_2[unit_length == 3]
#' exstra_wgs_pcr_2[unit_length == 3]$db$locus
#' 
#' # Extract all coding repeats
#' exstra_wgs_pcr_2[gene_region == "coding"]
#' exstra_wgs_pcr_2[gene_region == "coding"]$db$locus
#' 
#' # Extract all case samples:
#' exstra_wgs_pcr_2[, group == "case"]
#' 
#' # Extract by index 
#' exstra_wgs_pcr_2[2, 5]
#' 
#' @export
`[.exstra_score` <- function(x, loc, samp) {
  assert("locus is not the key of x$data", key(x$data)[1] == "locus")
  assert("sample is not the key of x$samples", key(x$samples)[1] == "sample")
  assert("locus not the key of x$db", key(x$db)[1] == "locus")
  if(!missing(loc)) {
    x$db <- x$db[eval(substitute(loc)), nomatch=0]
    setkey(x$db, locus)
  }
  if(!missing(samp)) {
    x$samples <- x$samples[eval(substitute(samp)), nomatch=0]
    setkey(x$samples, sample)
  }
  x$data <- x$data[x$db$locus][sample %in% x$samples$sample]
  setkey(x$data, locus, sample)
  x
}

#' Plot ECDF curves from an exstra_score object
#' 
#' 
#' @param x exstra_score object.
#' @param loci character; Only plot at the loci specified/ 
#'        No effect if NULL. 
#' @param sample_col Specify the colours for the given samples as a named character vector.
#'        names() are sample names, while the values are characters specifying colours.
#'        If specified, other samples are black by default, but may have an alpha value
#'        less than 1 for transparency. 
#' @param refline logical; if TRUE, then vertical lines are drawn indicating the reference size. 
#'        This does not take into account mismatches in the STR, or false repeat hits.
#' @param ylab Label for the y-axis. 
#' @param verticals logical; if TRUE, vertical lines are drawn at steps. See \code{\link{plot.stepfun}}.
#' @param pch numeric; the plotting character size.
#' @param xlim numeric of length 2; the x-limits for the plot. 
#' @param ylim numeric of length 2; the y-limits for the plot.
#' @param alpha_control Transparency alpha value for the control samples. 
#' @param alpha_case Transparency alpha value for the case samples. No effect if NULL.
#' @param xlinked If not "all", limits X chromosome loci to only "male" or "female" samples.
#'        Can also filter on only samples with "known" or "missing" sex. 
#' @param xlinked.safe logical; if TRUE when xlinked is filtering on "male" or "female", 
#'        an error will be raised if any samples have an undefined sex. 
#'        If FALSE, samples with missing sex will also be filtered when xlinked is "male" or "female". 
#' @param x_upper_missing If xlim is not specified, then loci with no data will have this
#'        upper x-axis limit. 
#' @param main Plot title.
#' @param ... Further arguments to the \code{plot.ecdf} function.
#' 
#' @examples 
#' plot(exstra_wgs_pcr_2["HD"])
#' 
#' plot(exstra_wgs_pcr_2, "HD", sample_col = c("WGSrpt_10" = "red", "WGSrpt_12" = "blue"))
#' 
#' # X-linked disorder, show male samples only:
#' plot(exstra_wgs_pcr_2, "SBMA", xlinked = "male")
#' 
#' @importFrom stats ecdf
#' @importFrom graphics abline grid text
#' @export
plot.exstra_score <- function(x, loci = NULL, sample_col = NULL, 
  refline = TRUE, ylab="Fn(x)", verticals = TRUE,
  pch = 19, xlim, ylim = c(0,1), alpha_control = 0.5, alpha_case = NULL, 
  xlinked = "all", xlinked.safe = TRUE, x_upper_missing = 150, 
  main = NULL, ...) {
  # Plot ECDFs of rep score data
  # sample_col should be a named vector, sample names as the name and color as the value
  # refline: if TRUE, include reference
  # xlinked: For loci on X chromosome, "all" for all samples, "male" and "female" for only that sex
  rsc <- x # for clarity
  if(is.null(loci)) {
    strlocis <- loci(rsc)
  } else {
    strlocis <- loci
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
      if(is.null(main)) {
        main.title <- paste(loci_text_info(rsc, locus.name), "score ECDF")
      } else {
        main.title <- main 
      }
      if(xlinked != "all" && grepl("X", rsc$db[locus.name]$inheritance)) {
        plot_data <- filter_sex(rsc, xlinked, xlinked.safe)$data[locus == locus.name]
        if(is.null(main)) {
          main.title <- paste(main.title, paste0(xlinked, 's'))
        }
      } else {
        plot_data <- rsc$data[locus.name, nomatch = 0]
      }
      if(missing(xlim)) {
        if(plot_data[, .N]) {
          xlim_1 <- c(0, max(plot_data$mlength))
        } else {
          xlim_1 <- c(0, x_upper_missing)
        }
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
        abline(v = loci_normal_exp(rsc, locus.name), col = c("blue", "red"), lty = 3:4)
      }
      black_trans <- rgb(0, 0, 0, alpha = alpha_control)
      if(is.null(sample_col)) {
        sample_col = ifelse(rsc$samples$group == 'case', rgb(1, 0, 0, alpha_case), black_trans)
        names(sample_col) <- rsc$samples$sample
      } 
      if(!is.null(alpha_case)) {
        sample_col <- add_alpha_(sample_col, alpha_case)
      }
      for(samp in c(setdiff(unique(plot_data$sample), names(sample_col)), intersect(unique(plot_data$sample), names(sample_col)))) {
        plot(ecdf(plot_data[sample == samp, rep]), add = TRUE, col = replace(sample_col[samp], is.na(sample_col[samp]), black_trans), verticals = verticals,
          pch = pch, ...)
      }
      if(plot_data[, .N] == 0) {
        text(x_upper_missing / 2, 0.5, "No data", cex = 2.5)
      }
    }
  }
}



# TODO: easy renaming of samples


# Copy an exstra_score object.
# 
# Prevents changing both objects on changes by reference, that do not copy on write. 
# @param x exstra_score object to copy.
# @export
copy.exstra_score <- function(x) {
  x$data %<>% copy()
  x$db %<>% copy()
  x$samples %<>% copy()
  x$input_type %<>% copy()
  x
}


# Number of data points in exstra_score object
# @param x exstra_score object.
# 
length.exstra_score <- function(x) {
  x$data[, .N]
}
  
  
#' @export
`length<-.exstra_score` <- function(x, value) {
  stop("Cannot reassign length to exstra_score object.")
}

# Dimension of exstra_score object
# 
# @param x An exstra_score object.
# 
# @return Vector of length two, the number of loci and number of samples respectively.
# 
dim.exstra_score <- function(x) {
  c(exstra_wgs_pcr_2$db[, .N], exstra_wgs_pcr_2$samples[, .N])
}

#' Convert a compatable object to the exstra_score class
#' 
#' @param x An object with a class that inherits from exstra_score
#' @param copy Should the object be copied? Slower if TRUE, but in-place data.table
#'             operations will not change both objects. This option is here to remind 
#'             users that this is normally copied by reference. 
#' @export
as.exstra_score <- function(x, copy = FALSE) {
  #
  assert("x should inherit from class exstra_score.", is.exstra_score(x))
  if(copy) {
    x <- copy.exstra_score(x)
  }
  structure(list(data = x$data, db = x$db, input_type = x$input_type, samples = x$samples), class = c("exstra_score", "exstra_db"))
}


# Verify keys of exstra_score
verify.exstra_score <- function(x) {
  assert("Object should inherit from class exstra_score.", is.exstra_score(x))
  assert("Key of x$data should be 'locus', 'sample'", identical(key(x$data), c("locus", "sample")))
  assert("Key of x$samples should be 'sample'.", key(x$samples) == "sample")
  verify.exstra_db(exstra_wgs_pcr_2)
}
