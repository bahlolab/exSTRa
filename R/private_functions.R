qq.unif <- function(ps, threshold = NULL, fudge = 2e-16, ...) {
  ps_in <- ps
  ps <- -log10(as.vector(ps) + fudge)
  p_null <-  -log10(ppoints(length(ps)) + fudge)
  threshold_log <- -log10(threshold + fudge)
  qqplot(p_null, ps, 
    xlab = expression(Expected~~-log[10](italic(p) + e)), 
    ylab = expression(Observed~~-log[10](italic(p) + e)),
    #ylim = range(0, ps),
    ylim = range(0, -log10(fudge)), # not quite what it should be
    ...)
  ps.order <- order(ps)
  abline(a = 0, b = 1, col = "red")
  if(!is.null(threshold)) {
    abline(h = threshold_log , col = "blue")
  }
  if(!is.null(names(ps_in))) {
    text(p_null %>% rev, ps[ps.order], labels = ifelse(threshold_log <= ps[ps.order], names(ps_in[ps.order]), ""), 
      pos = 2.5)
  }
}

remove_below_quant <- function(loc_data, quant = 0.5) {
  loc_data[, ord := order(prop), by = sample] [, keepme := ord > max(ord) * quant, by = sample]
  loc_data <- loc_data[identity(keepme), ]
  loc_data[, c("ord","keepme") := NULL]
  loc_data
}

plot_many_str_score <- function(strscore, typename, plot_cols, loci = NULL, 
  color_only = NULL, plottypes = 1, dirbase = "images/", 
  alpha_control = 0.5, alpha_case = NULL,
  legend = TRUE, legend_control = TRUE, controls_label = "controls", custom_legend = NULL,
  ...) {
  # typename a non-empty character string
  # plottypes should be vector of 1 to 3 can be 1:3
  # plot_cols may be a list, named by loci
  # color_only is a list, each item indicating the items to be coloured
  #TODO: check that sample names are correct
  # legend_custom, a named vector of colors for the legend
  if(any(is.element(plottypes, 2:3))) {
    dir.create(paste0(dirbase, typename), recursive = TRUE)
  }
  if(is.null(loci)) {
    loci <- loci(strscore)
  }
  for(i in plottypes) {
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
      plot(strscore, locus = loc, sample_col = plot_cols_this, 
        alpha_control = alpha_control, alpha_case = alpha_case, ...)
      plot_cols_this <- plot_cols_this[names(plot_cols_this) %in% strscore$samples$sample] # don't add legend for samples not shown
      leg_labels <- names(plot_cols_this)
      #if(length(plot_cols_this) != length(plot_cols)) {
      if(!is.null(alpha_case)) {
        plot_cols_this <- exSTRa::add.alpha(plot_cols_this, alpha_case)
      }
      if(legend_control && length(plot_cols_this) != strscore[loc]$samples[, .N]) {
        leg_labels <- c(leg_labels , controls_label)
        plot_cols_this <- c(plot_cols_this, rgb(0, 0, 0, alpha_control))
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

str_score_wmw <- function(strscore, quant = 0.5, B = NULL, loci = NULL, ncpus = 4, 
  t.test = FALSE, seed = 8109587) {
  # B = NULL for no permutation testing
  if(is.null(loci)) {
    loci <- loci(strscore)
  }
  set.seed(seed)
  score_ps <- matrix(nrow = length(strscore$samples$sample), ncol = length(loci))
  colnames(score_ps) <- loci
  rownames(score_ps) <- strscore$samples$sample
  for(loc in loci) {
    for(samp in strscore$samples$sample) {
      loc_data <- strscore$data[loc]
      loc_data[, highprop := quantile(prop, quant, names = FALSE), by = list(sample)]
      loc_data <- loc_data[prop > highprop]
      loc_data$this_group <- factor("control", levels = c("case", "control"))
      loc_data[sample == samp, this_group := "case"]
      if(t.test) {
        wt <- tryCatch(
          t.test(loc_data[this_group == "case", prop], 
            loc_data[this_group == "control", prop], 
            alternative = "greater")$p.value,
          error = function(x) { NA }
        )
      } else if(is.null(B)) {
        wt <- tryCatch(wilcox.test(loc_data[this_group == "case", prop], 
          loc_data[this_group == "control", prop], 
          alternative = "greater")$p.value,
          error = function(x) { NA }  
        )
      } else {
        wt <- tryCatch(wilcox_test(prop ~ this_group, 
          data = data.frame(loc_data), 
          alternative = "greater",
          distribution = approximate(B = B, ncpus = ncpus,  parallel = "multicore")
        ) %>% pvalue, # NOTE: this is possibly a bug, but only for perumation testing of Mann-Whitney
          error = function(x) { NA }
        )
      }
      
      score_ps[samp, loc] <- wt
    }
  }
  score_ps
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
  if(class(x) != "matrix") {
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

analyse_str_score_mw <- function(strscore, plot_cols, filebase = c(), actual = NULL, B = NULL, quant = 0.5, 
  loci = NULL, ncpus = 4, significance.threshold = NULL, 
  t.test = FALSE, cex.main.hist = 0.9, cex.main.qq = cex.main.hist) {
  # analyse the STR score with a Mann-Whitney or t.test
  if(t.test) {
    filebase <- c(filebase, "ttest")
  } else {
    filebase <- c(filebase, "mw")
    if(!is.null(B)) {
      filebase <- c(filebase, paste0("perm", B))
    }
  }
  if(quant != 0.5) {
    filebase <- c(filebase, paste0("q", quant))
  }
  if(!is.null(significance.threshold)) {
    filebase <- c(filebase, paste0("sigthres", signif(significance.threshold, 3)))
    bf_name <- ""
  } else {
    bf_name <- "-0.05bf"
  }
  filebase <- paste(filebase, collapse = "-")
  dirbase <- paste0("output/STR-score-", filebase, "/")
  dir.create(dirbase, recursive = TRUE)
  # Do everything we want to do, arrrrrr
  score_p <- str_score_wmw(strscore, B = B, quant = quant, loci = loci, ncpus = ncpus,
    t.test = t.test)
  
  # CSV files
  write.csv(score_p, paste0(dirbase, "score-comparison-", filebase, "-p.csv"))
  score_ps_summary <- summary_ps(score_p, significance.threshold = significance.threshold)
  write.csv(rbind(score_ps_summary$by.locus, 
    "Total" = list("Total", score_ps_summary$by.locus$n.sig %>% sum, NA)), 
    paste0(dirbase, "score-comparison-", filebase, "-summary-locus.csv"), row.names = FALSE)
  write.csv(rbind(score_ps_summary$by.sample, 
    "Total" = list("Total", score_ps_summary$by.sample$n.sig %>% sum, NA)), 
    paste0(dirbase, "score-comparison-", filebase, "-summary-sample.csv"), row.names = FALSE)
  
  # ECDF plots
  #TODO: fix up for significance.threshold
  plot_many_str_score(strscore, paste0("STR-score-", filebase, bf_name), plot_cols, 
    color_only = str_significant(score_p, significance.threshold), plottypes = 3, dirbase = dirbase)
  
  # Histograms
  dirbase_ext <- paste0(dirbase, "STR-score-", filebase, "-histograms/")
  dir.create(dirbase_ext, recursive = TRUE)
  pdf(paste0(dirbase_ext, "/all-histogram-", filebase, ".pdf"), useDingbats=FALSE)
  score_p %>% as.vector %>% hist(20, main = "Histogram of all loci")
  dev.off()
  
  for(loc in loci(strscore)) {
    pdf(paste0(dirbase_ext, "/histogram-", filebase, "-", loc, ".pdf"), useDingbats=FALSE)
    score_p[, loc] %>% as.vector %>% hist(20, 
      main = paste("Histogram of", loci_text_info(strscore, loc)),
      cex.main = cex.main.hist)
    #abline(v = 0.05 / length(score_ps), col = "blue")
    dev.off()
  }
  
  # Q-Q plots
  dirbase_ext <- paste0(dirbase, "/STR-score-", filebase, "-QQplots/")
  dir.create(dirbase_ext, recursive = TRUE)
  fudge <- ifelse(is.null(B), 2e-16, 0.5 / B)
  significance.threshold <- ifelse(is.null(significance.threshold), 0.05 / length(score_p), significance.threshold)
  pdf(paste0(dirbase_ext, "all-QQplots-", filebase, ".pdf"), useDingbats=FALSE)
  qq.unif(score_p, threshold = significance.threshold, main = "All loci", 
    fudge = fudge)
  dev.off()
  
  for(loc in colnames(score_p)) {
    pdf(paste0(dirbase_ext, "QQplot-", filebase, "-", loc, ".pdf"), useDingbats=FALSE)
    qq.unif(score_p[, loc], threshold = significance.threshold, 
      main = paste("Q-Q plot", strloci_text_info(strscore, loc)), fudge = fudge, 
      cex.main = cex.main.qq)
    dev.off()
  }
  
  # tidy it
  score_p_tidy <- score_p %>% as.data.frame()
  score_p_tidy$sample <- rownames(score_p_tidy)
  score_p_tidy %<>% gather(locus, pvalue, -sample)
  score_p_tidy <- data.table(score_p_tidy)
  score_p_tidy$expansion <- FALSE
  
  # Set the correct sample loci as expanded
  if(!is.null(actual)) {
    assert("actual is not a data.table (or NULL)", is.data.table(actual))
    setkey(score_p_tidy, sample, locus)
    setkey(actual, sample, locus)
    score_p_tidy[actual, expansion := TRUE]
    
    # Output ROC curves
    
    # pROC, not used 
    # roc(expansion ~ pvalue, score_p_tidy, plot = TRUE, ci = TRUE, colorize = TRUE)
    
    # ROCR
    preds <- with(score_p_tidy, prediction( predictions = -log10(pvalue), labels = expansion ))
    perf <- performance(preds, measure="tpr", x.measure="fpr")
    auc_value <- auc(expansion ~ pvalue, score_p_tidy) %>% as.numeric()
    
    pdf(paste0(dirbase, "/ROC-", filebase, ".pdf"), useDingbats=FALSE)
    plot(perf, colorize=T, lwd=2.5, main = paste0("ROC-", filebase))
    text(0.8, .15, labels = paste("AUC =", round(auc_value, 4)))
    dev.off()
    # performance(preds, measure="auc")
    
    # auc for pROC
    
  } else {
    auc_value <- NULL
  }
  
  list(score_p = score_p, score_p_tidy = score_p_tidy, auc = auc_value, dirbase = dirbase)
}


analyse_str_score_ps <- function(...) {
  # give the p values of analysing the STR score by different statistical tests
  # recommend that actual is left to NULL
  analyse_str_score_mw(..., actual = NULL)$score_p
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
  # method "al1" uses Munoz Rueda's AL1 algorithm to impute up to the desired number of data points
  #       If the number of data points is less than avaliable for a sample, then it is instead downsampled
  # method "quantile8" uses the quantile method type 8 (or number inferred)
  # method "al1_all" uses Munoz Rueda's AL1 algorithm to impute all data points to give quantiles
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
  assert("Cannot set both probs and n.quantiles", is.null(probs) || is.null(n.quantiles))
  if(is.null(probs) && is.null(n.quantiles)) {
    n.quantiles <- loc_data[, .N, by = sample][, quantile(N, read_count_quant, names = FALSE)] %>% round
  }
  method <- tolower(method)
  if(method == "al1" || method == "al1_all") {
    if(!is.null(probs)) {
      stop("probs cannot be manually set with al1 or al1_all")
    }
    probs <- seq(1 / n.quantiles, 1, length.out = n.quantiles)
  } else {
    if(is.null(probs)) {
      #probs <- ppoints(n.quantiles, 1/2)
      probs <- seq(0, 1, length.out = n.quantiles)
    } else {
      n.quantiles <- length(probs) # replace n.quantiles
    }
  }
  if(method == "quantile") {
    stop('"Please choose type for quantile with method = "quantile#"')
  }
  if(grepl("^quantile\\d", method)) {
    quantile_type <- sub("quantile", "", method) %>% as.numeric()
    method <- "quantile"
  }
  assert("samples is not the key of strscore$samples", key(strscore$samples)[1] == "sample")
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
    } else if(method == "al1") {
      if(n.quantiles < length(y)) {
        # sample down
        v <- sample_safe(y, n.quantiles)
      } else if (n.quantiles > length(y)) {
        # impute up
        v <- munoz_rueda_al1_include(y, n.quantiles - length(y))
      } else {
        # exact match
        v <- y
      }
      v <- sort(v)
    } else if(method == "al1_all") {
      v <- munoz_rueda_al1(y, n.quantiles) %>% sort
    } else {
      stop("Undefined method ", method)
    }
    quant.matrix[sampi, ] <- v
  }
  rownames(quant.matrix) <- sample
  # remove NA rows
  low.count.samples <- apply(quant.matrix, 1, function(x) { is.na(x) %>% all })
  quant.matrix <- quant.matrix[!low.count.samples, ]
  list(x = probs, y.mat = quant.matrix, low.count = names(low.count.samples[low.count.samples]))
}

trim_vector <- function(dim1, trim) {
  # takes an interger and trim value, gives indicies to keep
  # dim1: integer
  # trim: trimming proportion
  seq(ceiling(dim1 * trim) + 1, floor(dim1 * (1 - trim)), 1)
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
  assert("qmmat is not a matrix or list with $y.mat", is.matrix(qmmat))
  assert("sample is not a character or numeric", is.character(sample) || is.numeric(sample))
  assert("qs is not numeric", is.null(qs) || is.numeric(qs))
  assert("quant is not a single numeric from 0 and less than 1", is.numeric(quant), length(quant) == 1, quant >= 0, quant < 1)
  assert("qs_trim is not numeric of length 1", is.numeric(qs_trim), length(qs_trim) == 1)
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
  assert("trim is not numeric of length 1", is.numeric(trim), length(trim) == 1)
  assert("trim is not within 0 (inclusive) to 0.5", trim >= 0, trim < 0.5)
  if (trim > 0.3) {
    #TODO: make this a better test, such as the number of samples left (maybe I've already done this?) - Rick
    warning("trim is set to ", trim, ", removing ", trim * 200, "% of the data.")
  }
  assert("Maximum qs to output must be less or equal to than the number of xs - qs_trim", max(qs) <= dim(qmmat)[2] - qs_trim)
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
  assert("Trimming of extreme control samples set too high to keep at least two samples for variance estimation", length(trimseq) >= 2) 
  
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

quant_statistic_sampp <- function(qmmat, sample = NULL, qs = NULL, ...) {
  # get the quantile statistic for multiple samples
  # qmmat: a quantile matrix from the make_quantiles_matrix() function
  # sample: samples to get the statistic of. If NULL, give all samples in qmmat
  # qs: keeps the top number of quantiles, you probably do not want to use this
  # ... further arguments to quant_statistic(), most interesting is:
  #     quant keeps quantiles above its value only when qs is not specified. The default 
  #         keeps all values above the median (quant = 0.5)
  #     trim removes the top proportion of data points at each quantile 
  #         (rounded to the nearest sample) (0 <= trim < 1), default 0 at time of writing
  if(is.list(qmmat) && !is.null(qmmat$y.mat)) {
    qmmat <- qmmat$y.mat
  }
  assert("qmmat is not a matrix or list with $y.mat", is.matrix(qmmat))
  assert("sample is not a character, numeric or null", is.null(sample) || is.character(sample) || is.numeric(sample))
  assert("qs is not numeric", is.null(qs) || is.numeric(qs))
  assert("qs is not single", is.null(qs) || length(qs) == 1)
  if(is.null(sample)) {
    sample <- seq_len(dim(qmmat)[1])
  }
  t_out <- rep(NA, length(sample))
  ti <- 0
  for(samp in sample) {
    ti <- ti + 1
    t_out[ti] <- quant_statistic(qmmat, sample = samp, qs = qs, ...)
  }
  names(t_out) <- rownames(qmmat)
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
  assert("Trimming of extreme control samples set too high to keep at least two samples for variance estimation", length(trimseq) >= 2) 
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
        "simulate_quantile_matrix", 
        "quant_statistic",
        "trim_vector",
        "munoz_rueda_al1",
        "munoz_rueda_al1_include" 
      ),
      envir = environment()
    )
    # Run replicate
    sim_T <- parReplicate(cluster, B, simulate_quant_statistic())
    if(cluster_stop) {
      stopCluster(cluster)
    }
  } else {
    # not performing in parallel
    sim_T <- rep(0, B) # Don't let p-value get to zero
    for(i in seq_len(B)) {
      sim_T[i] <- simulate_quant_statistic()
    }
  }
  
  
  
  comp_p_function <- ecdf(sim_T)
  
  # p-values
  T_stat <- quant_statistic_sampp(qmmat, trim = trim, subject_in_background = subject_in_background, use_truncated = use_truncated_sd_in_quant_statistic) # TODO add ...
  # p <- 1 - comp_p_function(T_stat) + 1 / (B + 1)
  p <- rep(1.1, N)
  # calculate p-values
  for(i in seq_len(N)) {
    p[i] = c(sum(sim_T > T_stat[i]) + 1) / (B + 1)
  }
  names(p) <- names(T_stat)
  
  # generate ECDF and output
  list(ecdf = comp_p_function, p = p, sim_T = sim_T, T_stat = T_stat)
}



performance_from_stat_matrix <- function(stat_matrix, actual, main = NULL, pdfout = NULL,
  plot = TRUE, transform = identity) {
  # Generating performance from a statistics matrix, where
  # the rows are over samples, and columns loci. 
  # The known status of samples must also be provided
  #
  # stat_matrix: a matrix as specified above. 
  # actual: a data.table with columns sample and locus, consisting of conditional positives only
  # main: title for the ROC curve
  # pdfout, if specified, the plot is written to this path in PDF format
  # transform: a function to apply to the stat values
  
  # TODO: allow data.table() input, with sample and locus, as this will make this better
  
  # convert input matrix to a data.table()
  stat_dt <- stat_matrix %<>% reshape2:::melt.matrix(value.name = "stat", varnames = c("sample", "locus")) %>% data.table()
  
  # set expansion condition
  stat_dt$expansion_condition <- FALSE
  setkey(stat_dt, sample, locus)
  setkey(actual, sample, locus)
  stat_dt[actual, expansion_condition := TRUE] 
  
  # ROC curve
  preds <- with(stat_dt, prediction( predictions = transform(stat), labels = expansion_condition ))
  perf <- performance(preds, measure="tpr", x.measure="fpr")
  
  auc_value <- auc(expansion_condition ~ stat, stat_dt) %>% as.numeric()
  
  if(plot) {
    if(!is.null(pdfout)) {
      pdf(pdfout, useDingbats=FALSE)
    }
    plot(perf, colorize=TRUE, lwd=2.5, main = main) #, main = paste0("ROC-", filebase))
    text(0.8, .15, labels = paste("AUC =", round(auc_value, 4)))
    if(!is.null(pdfout)) {
      dev.off()
    }
  }
  
  list(auc = auc_value, stat_dt = stat_dt)
}

# summarise sensitivity/specificity 
sen_spec_all <- function(ps_list, truelist) {
  sen_spec <- data.frame(test = character(), sensitivity = numeric(), specificity = numeric(), tp = numeric(), fp = numeric(), fn = numeric(), tn = numeric(), stringsAsFactors = FALSE)
  for(px_name in names(ps_list)) {
    px <- ps_list[[px_name]]
    ss <- sensitivity_specificity(px, truelist)
    sen_spec <- rbind(sen_spec, list(test = px_name, sensitivity = ss$sensitivity, specificity = ss$specificity, tp = ss$confusion[1,1], fp = ss$confusion[2,1], fn = ss$confusion[1,2], tn = ss$confusion[2,2]))
    sen_spec$test %<>% as.character()
    #sen_spec <- rbind(sen_spec, c(sensitivity=ss$sensitivity, ss$specificity, (ss$confusion %>% as.vector)))
  }
  sen_spec
}

# 
midpoint_removal_imputing <- function(strscore, method, sort.in.original = TRUE, boxplot = FALSE) {
  #
  # sort.in.original: if TRUE, put the imputed values together with the orignal and sort
  # before calculating differences
  assert("Not strdata", is.strdata(strscore))
  quantile_diff <- data.table()
  for(loc in loci(strscore)) {
    for(samp in strscore$samples$sample) {
      loc_data <- strscore[loc, samp]
      if(loc_data$data[, .N] < 5) {
        # Since we need at least 3 points remaining for make_quantiles_matrix() to give an answer
        # 5 - 2 = 3 
        next
      }
      if(mod(loc_data$data[, .N], 2) == 0) {
        # even number of rows so remove lowest one
        loc_data$data <- loc_data$data[2:.N]
        #TODO: if 0 rows, this is bad
      }
      prune_data <- loc_data
      prune_data$data <- prune_data$data[order(rep)][seq(1, .N, 2)]
      n_out <- prune_data$data[, .N] - 1
      if(grepl("al1", method)) {
        qm_mid <- make_quantiles_matrix(prune_data, method = method, n.quantiles = n_out)
      } else {
        probs <- seq(1, 2 * n_out, 2) / (2 * n_out)
        #seq(3 / (n_out * 4), (n_out * 4 - 3) / (n_out * 4), length.out = n_out)
        qm_mid <- make_quantiles_matrix(prune_data, method = method, probs = probs)
      }
      
      # Consider that the sorting order may need to be different
      if(length(qm_mid$low.count) != 0) {
        diffs <- c()
      } else if(sort.in.original){
        newy <- sort(c(qm_mid$y.mat, prune_data$data$rep))
        diffs <- as.vector(newy - loc_data$data[order(rep)]$rep)
        diffs <- diffs[seq(2, length(diffs), 2)]
      } else {
        diffs <- as.vector(qm_mid$y.mat - loc_data$data[order(rep)][seq(2, .N, 2)]$rep)               
      }
      
      if(length(diffs) > 0) {
        if(length(loc_data$data[order(rep)][seq(2, .N, 2)]$rep) != length(qm_mid$y.mat)) {
          show(loc)
          show(samp)
          stop("About to insert things that are a different length")
        }
        quantile_diff <- rbind(quantile_diff, data.table(locus = loc, sample = samp, diff = diffs, qmid = as.vector(qm_mid$y.mat), probs = qm_mid$x))
      }
    }
    if(boxplot){
      boxplot(diff ~ sample, quantile_diff[locus == loc], main = strloci_text_info(strscore, loc))
      abline(h = 0, col = rgb(216, 98, 18, 150, maxColorValue = 255), lty = 4)
    }
  }
  
  quantile_diff$locus <- ordered(quantile_diff$locus, unique(quantile_diff$locus))
  quantile_diff
}



variance_of_trimmed <- function(x, fudge = 1) {
  # estimates the original variance of a trimmed observation
  # Does NOT perform trimming
  # Using Wikipedia variance (TODO: verify)
  # https://en.wikipedia.org/w/index.php?title=Truncated_normal_distribution&oldid=779765078
  # variance of trimmed  = <math>\sigma^2\left[1+\frac{\alpha\phi(\alpha)-\beta\phi(\beta)}{Z}
  # where <math>Z=\Phi(\beta)-\Phi(\alpha)</math>
  # \sigma^2 is underlying normal dist variance
  #
  # We want to reverse this!
  # Using same notation as Wikipedia, though the use of Z is kind of confusing
  # Note that also, this is for a truncated normal, not a trimmed where we fix the number of samples to remove
  # This is why we have a fudge parameter, that decreases alpha/increases beta a little 
  # fudge = 0 means no change, fudge = 1 means the change is 1/n * difference of alpha and beta
  
  # Centre data on 0, to work with formula
  y <- x - mean(x)
  
  if(var(y) == 0) {
    # No variance of input data, so return 0 to save numerical errors        
    return(0)
  }
  
  # These are probably incorrect as the actual alpha and beta should be
  alpha <- min(y)
  beta <- max(y) # this
  ab_diff <- beta - alpha
  
  # fudging
  alpha <- alpha - ab_diff * fudge / length(x)
  beta <- beta + ab_diff * fudge / length(x)
  
  # Preprepare variables
  Z <- pnorm(beta) - pnorm(alpha)
  phi_alpha <- dnorm(alpha)
  phi_beta <- dnorm(beta)
  
  # output 
  vtrim <- var(y) / 
    (
      1 +
        (alpha * phi_alpha - beta * phi_beta) / Z +
        ( (phi_alpha - phi_beta ) / Z ) ^ 2
    )
  if(is.nan(vtrim)) {
    # We divided by zero as there was no data variance
    vtrim <- 0
  }
  if(vtrim < 0) {
    # rounding error may make this happen
    if(vtrim <- 0.01) {
      stop("Potential exSTRa software bug, vtrim less than 0.01, vtrim = ", vtrim)
    }
    vtrim <- 0
  }
  vtrim
}

sd_of_trimmed <- function(...) {
  # Convienience function for variance_of_trimmed()
  # gives the square root of variance_of_trimmed()
  sqrt( variance_of_trimmed(...) )
}


# This runs the simulation with the given paramets
# Use ... to pass options to simulate_ecdf_quant_statistic()
# This also passes remaining options onto quant_statistic() and quant_statistic_sampp()
simulation_run <- function(data, B = 99, trim = 0.15, ...) {
  N <- data$samples[, .N]
  L <- data$db[, .N]
  p.matrix <- matrix(nrow = N, ncol = L)
  rownames(p.matrix) <- data$samples$sample
  colnames(p.matrix) <- loci(data)
  qmmats <- list()
  xecs <- list()
  for(loc in colnames(p.matrix)) {
    cat("Simulating distribution for", loc, "\n")
    qm_loop <- make_quantiles_matrix(data, loc, read_count_quant = 1, 
      method = "quantile7", min.n = 3)
    xec <- simulate_ecdf_quant_statistic(qm_loop, B, ...)
    p.matrix[names(xec$p), loc] <- xec$p
    qmmats[[loc]] <- qm_loop
    xecs[[loc]] <- xec
  }
  
  list(p.matrix = p.matrix, qmmats = qmmats, xecs = xecs)
}

# parallel replicate, from:
#https://stackoverflow.com/questions/19281010/simplest-way-to-do-parallel-replicate
# usage: 
# cl <- makePSOCKcluster(3) # set 3 to size of cluster (number of cores)
# hist(parReplicate(cl, 100, mean(rexp(10))))
parReplicate <- function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE)
  parSapply(cl, integer(n), function(i, ex) eval(ex, envir=.GlobalEnv),
    substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)

