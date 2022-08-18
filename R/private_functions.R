remove_below_quant <- function(loc_data, quant = 0.5) {
  loc_data[, ord := order(prop), by = sample] [, keepme := ord > max(ord) * quant, by = sample]
  loc_data <- loc_data[identity(keepme), ]
  loc_data[, c("ord","keepme") := NULL]
  loc_data
}


#' @importFrom grDevices dev.off pdf
plot_many_str_score <- function(strscore, typename, plot_cols, loci = NULL, 
  color_only = NULL, plot_types = 1, dirbase = "images/", 
  alpha_control = 0.5, alpha_case = NULL,
  legend = TRUE, legend_control = TRUE, controls_label = "controls", 
  cases_label = NULL, custom_legend = NULL,
  ...) {
  # typename a non-empty character string
  # plot_types should be vector of 1 to 3 can be 1:3
  # plot_cols may be a list, named by loci
  # color_only is a list, each item indicating the items to be coloured
  #TODO: check that sample names are correct
  # legend_custom, a named vector of colors for the legend
  if(any(is.element(plot_types, 2:3))) {
    # dir.create(paste0(dirbase, typename), recursive = TRUE)
    # The version below is more intuitive for other users
    dir.create(dirbase, recursive = TRUE)
    
  }
  if(is.null(loci)) {
    loci <- loci(strscore)
  }
  for(i in plot_types) {
    if(i == 2) { pdf(paste0(dirbase, typename, ".pdf"), useDingbats=FALSE) }
    for(loc in loci) {
      if(i == 3) { 
        pdf(paste0(dirbase, typename, "/", loc, "-", typename,".pdf"), useDingbats=FALSE) 
        par(mar = c(4, 4, 3, 0) + 0.1)
      }
      if(is.list(plot_cols)) {
        plot_cols_this <- plot_cols[[loc]]
      } else {
        plot_cols_this <- plot_cols
      }
      if(!is.null(color_only)) {
        if(is.list(color_only)) {
          plot_cols_this <- plot_cols_this[color_only[[loc]]]
        } else {
          stop("color_only should be list")
        }
      }
      plot(strscore, loci = loc, sample_col = plot_cols_this, 
        alpha_control = alpha_control, alpha_case = alpha_case, ...)
      plot_cols_this <- plot_cols_this[names(plot_cols_this) %in% strscore$samples$sample] # don't add legend for samples not shown
      leg_labels <- names(plot_cols_this)
      #if(length(plot_cols_this) != length(plot_cols)) {
      if(!is.null(alpha_case)) {
        plot_cols_this <- add_alpha_(plot_cols_this, alpha_case)
      }
      if(legend_control && length(plot_cols_this) != strscore[loc]$samples[, .N]) {
        leg_labels <- c(leg_labels, controls_label)
        plot_cols_this <- c(plot_cols_this, rgb(0, 0, 0, alpha_control))
      }
      if(!is.null(cases_label)) {
        leg_labels <- c(leg_labels, cases_label)
        plot_cols_this <- c(plot_cols_this, rgb(1, 0, 0, alpha_case))
      }
      if(legend) {
        if(is.null(custom_legend)) {
          legend("bottomright", leg_labels, lty = 1, pch = 16, col = plot_cols_this, bg = "white", cex = 0.77)
        } else {
          legend("bottomright", names(custom_legend), lty = 1, pch = 16, col = custom_legend, bg = "white", cex = 0.77)
        }
      }
      if(i == 3) { dev.off() }
    }
    if(i == 2) { dev.off() }
  }
}


sample_safe <- function(x, size, replace = FALSE, ...) {
  # Sometimes the input is only of length one, causing different behaviour
  if(length(x) == 1) {
    if(missing(size)){
      return(x)
    } else {
      if(replace) {
        return(rep(x, size))
      } else {
        stop("cannot take a sample larger than the population when 'replace = FALSE'")
      }
    }
  }
  sample(x, size, replace = replace, ...)
}

#' @importFrom stats quantile
make_quantiles_matrix <- function(strscore, loc = TRUE, sample = NULL, read_count_quant = 1, 
  probs = NULL, quant = NULL, method = "midquantile", n.quantiles = NULL, min.n = 3) {
  # Create the matrix of quantiles
  # read_count_quant: quantile of the read count for chosing the number of quantiles in the output
  # probs, quantile data points, only for method "midquantile" and is prioritised over other methods to choose quantile values
  # quant, remove the data points below the given quantile for each sample at each locus, this probably should not be used, instead derived from the statistic
  #
  # method "midquantile" uses Qtools to find the mid-quantile at ppoints(n, a=1/2)
  # method "quantile8" uses the quantile method type 8 (or number inferred)
  # n.quantiles sets the number of quantiles in the output matrix
  # min.n is the minimum number of observations for a sample at a locus to go into the matrix
  
  setkey(strscore$data, locus, sample)
  loc_data <- strscore$data[loc]
  if(identical(TRUE, loc) && loc_data[, .N, by = locus][, .N] != 1) {
    warning("More than one locus has been passed to make_quantiles_matrix with default loc value of TRUE")
  } 
  if(loc_data[, .N] == 0) {
    stop(paste("No data left after filtering locus"))
  }
  if(!is.null(quant)) {
    loc_data <- remove_below_quant(loc_data, quant)
    warning("Setting the quantile in make_quantiles_matrix() is not recommended")
  }
  if(is.null(sample)) {
    sample <- strscore$samples$sample
  }
  testit::assert("Cannot set both probs and n.quantiles", is.null(probs) || is.null(n.quantiles))
  if(is.null(probs) && is.null(n.quantiles)) {
    n.quantiles <- round(loc_data[, .N, by = sample][, quantile(N, read_count_quant, names = FALSE)])
  }
  method <- tolower(method)
    if(is.null(probs)) {
      #probs <- ppoints(n.quantiles, 1/2)
      probs <- seq(0, 1, length.out = n.quantiles)
    } else {
      n.quantiles <- length(probs) # replace n.quantiles
    }
  if(method == "quantile") {
    stop('"Please choose type for quantile with method = "quantile#"')
  }
  if(grepl("^quantile\\d", method)) {
    quantile_type <- as.numeric(sub("quantile", "", method))
    method <- "quantile"
  }
  testit::assert("samples is not the key of strscore$samples", key(strscore$samples)[1] == "sample")
  quant.matrix <- matrix(numeric(), length(sample), length(probs))
  for(sampi in seq_along(sample)) {
    # note this "sample" is the one within the object
    y <- loc_data[sample == strscore$samples[sampi, sample]]$rep
    if(length(y) < min.n) {
      # sometimes there is no data to impute
      v <- rep(as.numeric(NA), n.quantiles)
    } else if(method == "midquantile") {
      if(length(unique(y)) == 1) {
        # midquantile is no good if all the values are the same
        # TODO: this can probably be improved, can result in all quantiles being 
        #       the same as they should be for quantiles, but may be poor for 
        #       analysis
        v <- rep(y[1], n.quantiles)
      } else {
        xmid <- midquantile(y, probs = probs)
        v <- xmid$y
      }
    } else if(method == "quantile") {
      v <- quantile(y, probs, names = FALSE, type = quantile_type)
    } else {
      stop("Undefined method ", method)
    }
    quant.matrix[sampi, ] <- v
  }
  rownames(quant.matrix) <- sample
  # remove NA rows
  low.count.samples <- apply(quant.matrix, 1, function(x) { all(is.na(x)) })
  quant.matrix <- quant.matrix[!low.count.samples, ]
  list(x = probs, y.mat = quant.matrix, low.count = names(low.count.samples[low.count.samples]))
}

trim_vector <- function(dim1, trim) {
  # takes an interger and trim value, gives indicies to keep
  # dim1: integer
  # trim: trimming proportion
  ti <- trim_index_(dim1, trim)
  seq(ti[1], ti[2], 1)
}

# gives indexes for trimming
trim_index_ <- function(dim1, trim) {
  c(ceiling(dim1 * trim) + 1, floor(dim1 * (1 - trim)))
}


# parallel replicate, from:
#https://stackoverflow.com/questions/19281010/simplest-way-to-do-parallel-replicate
# usage: 
# cl <- makePSOCKcluster(3) # set 3 to size of cluster (number of cores)
# hist(parReplicate(cl, 100, mean(rexp(10))))
parReplicate <- function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE)
  parSapply(cl, integer(n), function(i, ex) eval(ex, envir=.GlobalEnv),
    substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)

