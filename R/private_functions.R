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


str_significant <- function(ps, significance.threshold = NULL) {
  if(is.null(significance.threshold)) {
    significance.threshold <- 0.05 / length(ps %>% as.vector)
  }
  sig_list <- list()
  for(loc in colnames(ps)) {
    sig_list[[loc]] <- which(ps[, loc] <= significance.threshold) %>% names
  }
  sig_list
}

summary_ps <- function(ps, significance.threshold = NULL) {
  str_sig <- str_significant(ps, significance.threshold)
  nloc <- str_sig %>% length
  by.locus <- data.frame(
    locus = str_sig %>% names, 
    n.sig =  sapply(str_sig, length), 
    significant = sapply(str_sig, function(x) {paste(x, collapse = ",")}),
    stringsAsFactors = FALSE
  )
  by.sample <- data.frame(
    row.names = rownames(ps), 
    sample = rownames(ps),
    n.sig =  0, 
    significant = rep("", dim(ps)[1]),
    stringsAsFactors = FALSE
  )
  by.sample$significant <- as.character(by.sample$significant) # because it kept being made into a factor :(
  for(rowi in seq_along(str_sig)) {
    for(samp in str_sig[[rowi]]) {
      if(by.sample[samp, "significant"] == "") {
        by.sample[samp, "significant"] <- names(str_sig)[rowi]
      } else {
        by.sample[samp, "significant"] <- paste(by.sample[samp, "significant"], names(str_sig)[rowi], sep = ",")
      }
      by.sample[samp, "n.sig"] %<>% `+`(1)
    }
  }
  list(by.locus = by.locus, by.sample = by.sample)
}

sensitivity_specificity <- function(x, truelist, significance.threshold = NULL) {
  if(inherits(x, "matrix")) {
    stop("x is not matrix")
  }
  if(is.data.table(truelist)) {
    # convert to the expected list format
    truelist_dt <- truelist
    truelist <- list()
    for(loc in colnames(x)) {
      truelist[[loc]] <- actual_expansions[locus == loc & sample %in% rownames(x), sample]
    }
  }
  strsig <- str_significant(x, significance.threshold)
  x_n <- length(x)
  confusion <- matrix(0, ncol = 2, nrow = 2)
  # for(loc in names(truelist)) {
  for(loc in colnames(x)) {
    trues <- truelist[[loc]] %>% length
    true_p <- is.element(truelist[[loc]], strsig[[loc]]) %>% sum
    false_n <- trues - true_p
    false_p <- (!is.element(strsig[[loc]], truelist[[loc]])) %>% sum
    confusion %<>% `+`(matrix(c(true_p, false_p, false_n, 0), ncol = 2, nrow = 2))
  }
  confusion[2,2] <-  x_n - sum(confusion)
  sensitivity <- confusion[1,1] / sum(confusion[1,])
  specificity <- confusion[2,2] / sum(confusion[2,])
  list(sensitivity = sensitivity, specificity = specificity, confusion = confusion)
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

quant_statistic <- function(qmmat, sample = 1, quant = 0.5, trim = 0.15, 
  qs = NULL, qs_trim = 0, 
  test = t.test, 
  subject_in_background = TRUE, 
  use_truncated_sd = FALSE, # maybe should be true
  median_mad = FALSE, # use median and mad instead of mean and sd
  fudge = 1, # fudge factor for trimmed SD
  ...) {
  # sample may be an integer indicating the sample index, or the sample name
  # quant keeps quantiles above its value only when qs is not specified. The default 
  #       keeps all values above the median (quant = 0.5)
  # trim removes the top proportion of data points at each quantile 
  #      (rounded to the nearest sample) (0 <= trim < 0.5, though close to 0.5 often will not work)
  # qs keeps the top number of quantiles
  # qs_trim trims the top quantiles by this number of quantiles
  # test is the test function to be used
  if(is.list(qmmat) && !is.null(qmmat$y.mat)) {
    qmmat <- qmmat$y.mat
  }
  testit::assert("qmmat is not a matrix or list with $y.mat", is.matrix(qmmat))
  testit::assert("sample is not a character or numeric", is.character(sample) || is.numeric(sample))
  testit::assert("qs is not numeric", is.null(qs) || is.numeric(qs))
  testit::assert("quant is not a single numeric from 0 and less than 1", is.numeric(quant), length(quant) == 1, quant >= 0, quant < 1)
  testit::assert("qs_trim is not numeric of length 1", is.numeric(qs_trim), length(qs_trim) == 1)
  if(median_mad && use_truncated_sd) {
    stop("Options median_mad and use_truncated_sd are mutually exclusive.")
  }
  if(is.null(qs)) {
    qs <- dim(qmmat)[2] - qs_trim - floor(dim(qmmat)[2] * quant)
    if(qs < 1) {
      qs <- 1
      warning("qs set to 1, qs_trim may be too high for data or may be too few data points")
    }
  }
  testit::assert("trim is not numeric of length 1", is.numeric(trim), length(trim) == 1)
  testit::assert("trim is not within 0 (inclusive) to 0.5", trim >= 0, trim < 0.5)
  if (trim > 0.3) {
    #TODO: make this a better test, such as the number of samples left (maybe I've already done this?) - Rick
    warning("trim is set to ", trim, ", removing ", trim * 200, "% of the data.")
  }
  testit::assert("Maximum qs to output must be less or equal to than the number of xs - qs_trim", max(qs) <= dim(qmmat)[2] - qs_trim)
  #qmmat[sample, ]
  xindexes <- seq(dim(qmmat)[2] - max(qs) - qs_trim + 1, dim(qmmat)[2] - qs_trim, 1)
  ts <- rep(NA, length(xindexes))
  ti <- 0
  if(is.character(sample)) {
    sample <- which(rownames(qmmat) == sample)
  }
  
  bg_n <- dim(qmmat)[1] # number of background samples
  if(! subject_in_background) {
    # We are removing the subject sample, so trim only on background samples
    bg_n <- bg_n - 1
  } 
  trimseq <- trim_vector(bg_n, trim)
  testit::assert("Trimming of extreme control samples set too high to keep at least two samples for variance estimation", length(trimseq) >= 2) 
  
  for(xi in xindexes) {
    ti <- ti + 1
    #show(ti)
    # trim the most extreme values
    if(subject_in_background) {
      # Do not remove the subject for creating background observations (cx)
      cx <- qmmat[, xi] # no longer remove sample to keep same se, rely on trimming
    } else {
      # Remove the subject sample from the background
      cx <- qmmat[-sample, xi]
    }
    # the test statistic
    if(median_mad) {
      S <- mad(cx)
      if(S == 0) {
        t <- 0 # don't add to the T statistic
      } else {
        t0 <- (qmmat[sample, xi] - median(cx)) / S
      }
    } else {
      cx <- sort(cx)[trimseq] # trimming of background
      if(use_truncated_sd) {
        S  <- sd_of_trimmed(cx, fudge = fudge)
        t0 <- (qmmat[sample, xi] - mean(cx)) / S
      } else {
        t0 <- tryCatch(
          test(qmmat[sample, xi], 
            cx, 
            alternative = "greater", 
            var.equal = TRUE, ...)$statistic,
          error = function(x) { 0 } # TODO: check that having this as 0 instead of NA is ok
        ) 
      }
    }
    #show(t0)
    ts[ti] <- t0
  }
  #show(ts)
  cumsum(rev(ts))[qs] / qs # we divide by the number of t statistics being summed
}

quant_statistic_sampp <- function(qmmat, sample = NULL, qs = NULL, 
  case_samples = NULL,
  ...) {
  # get the quantile statistic for multiple samples
  # qmmat: a quantile matrix from the make_quantiles_matrix() function
  # sample: samples to get the statistic of. If NULL, give all samples in qmmat
  # qs: keeps the top number of quantiles, you probably do not want to use this
  # case_samples: if not NULL, only calculate for these samples in case-control setting.
  #               Other cases are excluded in each calculation.
  # ... further arguments to quant_statistic(), most interesting is:
  #     quant keeps quantiles above its value only when qs is not specified. The default 
  #         keeps all values above the median (quant = 0.5)
  #     trim removes the top proportion of data points at each quantile 
  #         (rounded to the nearest sample) (0 <= trim < 1), default 0 at time of writing
  if(is.list(qmmat) && !is.null(qmmat$y.mat)) {
    qmmat <- qmmat$y.mat
  }
  testit::assert("qmmat is not a matrix or list with $y.mat", is.matrix(qmmat))
  testit::assert("sample is not a character, numeric or null", is.null(sample) || is.character(sample) || is.numeric(sample))
  testit::assert("qs is not numeric", is.null(qs) || is.numeric(qs))
  testit::assert("qs is not single", is.null(qs) || length(qs) == 1)
  
  ti <- 0
  if(!is.null(case_samples)) {
    if(is.null(sample)) {
      sample <- case_samples
    }
    testit::assert("All case samples should be in qmmat", all(case_samples %in% rownames(qmmat)))
    t_out <- rep(NA, length(sample) - length(case_samples))
    control_samples <- rownames(qmmat)[! rownames(qmmat) %in% case_samples]
    for(samp in sample) {
      ti <- ti + 1
      # Get only the correct qmmat columns
      qmmat0 <- qmmat[c(samp, control_samples),]
      t_out[ti] <- quant_statistic(qmmat0, sample = samp, qs = qs, 
        subject_in_background = FALSE, ...) # TODO: maybe trim = 0?
    }
    names(t_out) <- sample
  } else {
    # No samples marked explicitly as cases
    if(is.null(sample)) {
      sample <- seq_len(dim(qmmat)[1])
    }
    t_out <- rep(NA, length(sample))
    for(samp in sample) {
      ti <- ti + 1
      t_out[ti] <- quant_statistic(qmmat, sample = samp, qs = qs, ...)
    }
    names(t_out) <- rownames(qmmat) # may be a bug if only some samples are required
  }
  t_out
}

generate_T_stat_performance <- function(strscore, actual, main = NULL, pdfout = NULL,
  quant = 0.5, trim = 0.2, plot = TRUE) {
  # Generating performance for the T stat
  # strscore the usual
  # actual a data.table with columns sample and locus, consisting of conditional positives only
  # main: title for the ROC curve
  # pdfout, if specified, the plot is written to this path in PDF format
  T_stats_list <- list()
  for(loc in loci(strscore)) {
    qm <- make_quantiles_matrix(strscore, loc = loc, sample = NULL, read_count_quant = 1, 
      method = "quantile7", min.n = 3)
    T_stats_loc <- quant_statistic_sampp(qm, quant = quant, trim = trim) # quant at default of 0.5
    T_stats_list[[loc]] <- data.table(sample = names(T_stats_loc), T = T_stats_loc)
  }
  T_stats <- rbindlist(T_stats_list, idcol = "locus")
  
  # set expansion condition
  T_stats$expansion_condition <- FALSE
  setkey(T_stats, sample, locus)
  setkey(actual, sample, locus)
  T_stats[actual, expansion_condition := TRUE] 
  
  # ROC curve
  preds <- with(T_stats, prediction( predictions = T, labels = expansion_condition ))
  perf <- performance(preds, measure="tpr", x.measure="fpr")
  
  auc_value <- auc(expansion_condition ~ T, T_stats) %>% as.numeric()
  
  if(plot) {
    if(!is.null(pdfout)) {
      pdf(pdfout, useDingbats=FALSE)
    }
    plot(perf, colorize=T, lwd=2.5, main = main) #, main = paste0("ROC-", filebase))
    text(0.8, .15, labels = paste("AUC =", round(auc_value, 4)))
    if(!is.null(pdfout)) {
      dev.off()
    }
  }
  
  list(auc = auc_value, T_stats = T_stats)
}



simulate_ecdf_quant_statistic <- function(qmmat, B = 9999, trim = 0.15, 
  subject_in_background = TRUE, 
  sort_sim_qm = TRUE, 
  use_truncated_sd = FALSE,
  truncated_sd_fudge = 1,
  use_truncated_sd_in_quant_statistic = FALSE,
  median_mad = FALSE, # use median and MAD of data for rnorm (instead of mean and sd)
  median_mad_in_quant_statistic = FALSE, 
  parallel = FALSE, # perform in parallel
  cluster = NULL, # cluster for parallel. If NULL, then make one
  cluster_n = 4, # cluster size if cl is NULL
  T_stat = NULL, # T_stat named vector precomputed
  ...) {
  # Derive a simulated ECDF, and return an "ecdf" class object
  # qmmat: a quantile matrix from the make_quantiles_matrix() function
  # B: the number of simulations, if i is an integer, setting B = 10^i - 1 means the 
  #     resulting p-value is limited to i decimal places without needing to round.
  # ... further arguments to quant_statistic(), most interesting is:
  #     quant keeps quantiles above its value only when qs is not specified. The default 
  #         keeps all values above the median (quant = 0.5)
  #     trim removes the top proportion of data points at each quantile 
  #         (rounded to the nearest sample) (0 <= trim < 1), default 0 at time of writing
  
  # Development options
  # subject_in_background: TRUE/FALSE, should the subject be included in the calculation of quant stat? #TODO change code here?
  # sort_sim_qm: should the simulated quantile matrix be sorted across quantiles by sample? Sampling by rnorm means it is no longer in order. 
  
  # Developer reminder:
  # quant_statistic <- function(qmmat, sample = 1, quant = 0.5, trim = 0.2, qs = NULL, qs_trim = 0, test = t.test, ...)
  # ... goes to t.test()
  
  ### The work ###
  ## trim ##
  # first, sort the rows
  
  qmmat.y.sort <- apply(qmmat$y.mat, 2, sort) # TODO maybe not working
  # so that we can trim the top and bottom:
  trimseq <- trim_vector(dim(qmmat.y.sort)[1], trim)
  testit::assert("Trimming of extreme control samples set too high to keep at least two samples for variance estimation", length(trimseq) >= 2) 
  if(length(trimseq) < 6) {
    # TODO: 6 chosen arbitrarily, maybe could be better
    warning("Trimming of outliers leaves only ", length(trimseq), " observations.")
  }
  qmmat.y.trim <- qmmat.y.sort[trimseq, ]
  
  # determine mu and se
  # devel: is it norm?
  # > qqnorm(qmmat.y.trim)
  # > qqnorm(qmmat.y.trim[, 55])
  # > qqnorm(qmmat.y.trim[, 60])
  # > qqnorm(qmmat.y.trim[, 65])
  
  
  N <- dim(qmmat.y.sort)[1] # number of samples
  M <- dim(qmmat.y.sort)[2] # number of quantiles
  mu_vec <- apply(qmmat.y.trim, 2, mean)
  
  if(median_mad && use_truncated_sd) {
    stop("median_mad and use_truncated_sd are mutually exclusive options")
  }
  
  if(median_mad) {
    # Calculate the mean and se from the median and MAD of the untruncated
    # data instead
    mu_vec <- apply(qmmat.y.sort, 2, median)
    se_vec <- apply(qmmat.y.sort, 2, mad) / qnorm(3/4)
  } else if(use_truncated_sd) {
    # Trying to find the original SD from the trimmed data
    sd_here <- function(x) { sd_of_trimmed(x, fudge = truncated_sd_fudge)}
    # stop dodgy numerical errors:
    qmmat.y.trim <- round(qmmat.y.trim, 10)
    se_vec <- apply(qmmat.y.trim, 2, sd_here)
  } else {
    se_vec <- apply(qmmat.y.trim, 2, sd)
  }
  
  # simulate B times
  # mu_vec and se_vec are already in this name space,
  # hence this function is only defined here
  sort_rev <- function(x) { 
    # NOT USED
    # Simple reversed sort
    sort(x, decreasing = TRUE)
  }
  simulate_quantile_matrix <- function() {
    sqm <- t (replicate(N, rnorm(M, mu_vec, se_vec)))
    # as this means the quantiles of a sample are no longer ordered, we (maybe) sort
    if(sort_sim_qm){
      for(i in seq_len(N)) {
        sqm[i, ] <- sort(sqm[i,  ])
      }
    }
    sqm
  }
  simulate_quant_statistic <- function() {
    simu <- simulate_quantile_matrix()
    do.call(quant_statistic, c(list(simu), triple_dots))
  }
  simulate_quant_statistic_sampp <- function() {
    simu <- simulate_quantile_matrix()
    do.call(quant_statistic_sampp, c(list(simu), triple_dots))
  }
  
  triple_dots <- c(list(...), list(subject_in_background = subject_in_background, 
    use_truncated_sd = use_truncated_sd_in_quant_statistic, 
    truncated_sd_fudge = truncated_sd_fudge,
    median_mad = median_mad_in_quant_statistic))
  
  # For some weird reason this didn't work:
  #sim_T <- replicate(B, simulate_quant_statistic())
  if(parallel) {
    if(is.null(cluster)) {
      cluster <- makePSOCKcluster(cluster_n)
      cluster_stop <- TRUE
    } else {
      cluster_stop <- FALSE
    }
    # Load required functions onto cluster nodes
    clusterEvalQ(cluster, { 
      library(exSTRa); 
      library(testit); 
      library(data.table);
      library(magrittr)
    })
    #put objects in place that might be needed for the code
    #... then parallel replicate...
    clusterExport(cluster,
      c("triple_dots",
        "subject_in_background", 
        "use_truncated_sd_in_quant_statistic", 
        "truncated_sd_fudge",
        "sort_sim_qm", 
        "N", "M", "mu_vec", "se_vec", 
        "simulate_quant_statistic",
        "simulate_quant_statistic_sampp",
        "simulate_quantile_matrix", 
        "quant_statistic",
        "trim_vector"
      ),
      envir = environment()
    )
    # Run replicate
    sim_T <- as.vector(
            parReplicate(cluster, B, simulate_quant_statistic_sampp())
    )
    if(cluster_stop) {
      stopCluster(cluster)
    }
  } else {
    # not performing in parallel
    # TODO here:
    sim_T <- rep(0, B*N) # Don't let p-value get to zero
    for(i in seq_len(B)) {
      sim_T[((i-1)*N + 1) : (i*N) ] <- simulate_quant_statistic_sampp()
      # TODO: insert stopping-code here
    }
  }
  
  
  # Following is not relevant to output
  # comp_p_function <- ecdf(sim_T)
  
  # p-values
  if(is.null(T_stat)) {
    T_stat <- quant_statistic_sampp(qmmat, trim = trim, subject_in_background = subject_in_background, use_truncated = use_truncated_sd_in_quant_statistic) # TODO add ...
  }
  # p <- 1 - comp_p_function(T_stat) + 1 / (B + 1)
  p <- rep(1.1, N)
  # calculate p-values
  for(i in seq_len(N)) {
    p[i] = (sum(sim_T > T_stat[i]) + 1) / (length(sim_T) + 1)
  }
  names(p) <- names(T_stat)
  
  # generate ECDF and output
  # list(ecdf = comp_p_function, p = p, sim_T = sim_T, T_stat = T_stat)
  list(p = p, sim_T = sim_T, T_stat = T_stat)
}


# parallel replicate, from:
#https://stackoverflow.com/questions/19281010/simplest-way-to-do-parallel-replicate
# usage: 
# cl <- makePSOCKcluster(3) # set 3 to size of cluster (number of cores)
# hist(parReplicate(cl, 100, mean(rexp(10))))
parReplicate <- function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE)
  parSapply(cl, integer(n), function(i, ex) eval(ex, envir=.GlobalEnv),
    substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)

