# The log linear tests

is.str_loglin_test <- function(x) inherits(x, "str_loglin_test")

#TODO: make this the actual new function we want
str_loglin_test_new <- function(data, loci, statistic, df, p.value, reads.total,  
                                    group_case,
                                    group_control,
                                    group_null)
{
  assert("Input should be data.frames", is.data.frame(statistic), is.data.frame(df), is.data.frame(p.value), is.data.frame(reads.total))
  assert("data should be of class strdata", is.strdata(data))
  assert("loci should be a character vector", is.character(loci))
  structure(
    list(data = data, loci = loci, statistic = statistic, df = df, p.value = p.value, reads.total = reads.total, 
         group_case = group_case,
         group_control = group_control,
         group_null = group_null), 
    class = c("str_loglin_test")
  )
}

strloci.str_loglin_test <- function(X) {
  strloci(X$data)
}

str_loglin_test <- function(data,
                                       cols, 
                                       keep.cols = cols,
                                       collapsing.guide = NULL, 
                                       loci = strloci(data), 
                                       require.nozero = FALSE,
                                       group_null = NULL,
                                       group_control = "control",
                                       group_case = "case",
                                       allow_uncounted_loci = FALSE,
                                       include.test.case = FALSE
) {
  # Runs a log-linear test on strcounts
  
  assert("data is required to be of the class strdata", inherits(data, "strdata"))
  assert("loci should be an atomic vector", is.atomic(loci))  
  
  # Pull the input data apart
  features <- data$data
  groups <- data$data$group
  # replace states by groups
  
  # Prepare inputs
  loci <- as.character(loci)
  
  # Warnings
  if(length(setdiff(loci, levels(features$locus))) > 0) {
    mess <- "loci are given that are not in the data."
    if(allow_uncounted_loci) {
      warning(mess)      
    } else {
      stop(mess, " Use option allow_uncounted_loci = TRUE if this is ok.")
    }
  }
  
  # Storing results
  disease.anova.chisq.statistics <- data.frame() #TODO: make data.table
  disease.anova.chisq.df <- data.frame()
  disease.anova.chisq.p.value <- data.frame()
  disease.read.total <- data.frame()
  for(locus.name in loci) {
    # Need to first sum the data over the normal samples and the 
    # expanded samples
    cat("STR: ", locus.name, "\n")
    
    # data
    locus.counts <- features[locus.name, c("sample", "group", keep.cols), with = F]
    setkey(locus.counts, sample)
        
    if(!is.null(collapsing.guide)) {
      stop("Collapsing not yet implemented") # TODO
      for(simple_feature in names(collapsing.guide)) {
        cont.table.1.sample <- collapse.contigency.table.single(cont.table.1.sample, 
                                                                grep(collapsing.guide[simple_feature], rownames(cont.table.1.sample), value = T), name = simple_feature)
      }
    }
    
    if(require.nozero) {
      #cont.table.1.sample <- cont.table.1.sample[apply(cont.table.1.sample, 1, sum) != 0, , drop = F] # Remove zero rows
      delete <- c()
      for(bin.i in seq(3, dim(locus.counts)[2], 1)) {
        if(is.null(group_null)) {
          total <- sum(locus.counts[group == group_control, bin.i, with = F])
        } else {
          total <- sum(locus.counts[group == group_null, bin.i, with = F])
        }
        if(total == 0) {
          delete <- c(delete, bin,i)
        }
      }
      locus.counts <- locus.counts[, -delete, with = F]
    }
    
    locus.counts.long <- reshape(locus.counts, direction = "long", varying = names(locus.counts)[-c(1,2)], v.names = "reads", times = names(locus.counts)[-c(1,2)], timevar = "bin")
    locus.counts.long$affected <- as.logical(NA)
    setkey(locus.counts.long, sample, bin)
    
    for(sample.name in data$samples$sample) {
      sample.data <- locus.counts[sample.name, nomatch = 0] # TODO: This was a problem, but I don't think it is anymore
      if(dim(sample.data)[1] == 0) {
        next
      }
      assert("Sample appears to be in multiple groups", length(unique(sample.data$group)) == 1)
      locus.counts.long$affected <- as.logical(NA)
      if(is.null(group_null)) {
        # No samples given for null distribution, so use the leave-one-out procedure for the null distribution
        locus.counts.long[group == group_control]$affected <- FALSE
        locus.counts.long[sample.name]$affected <- TRUE # turns out we do the same regardless of group
      } else {
        # Have samples for null distribution, so compare controls to the null and cases to the null
        # set: cont.table.1.sample.v
        stop("This feature not yet implemented for nulls, controls and cases")
      }
      if(include.test.case) {
        the.keys <- key(locus.counts.long)
        lc.case.pop <- locus.counts.long[sample.name]
        lc.case.pop$affected <- FALSE
        locus.counts.long <- rbind(locus.counts.long, lc.case.pop)
        setkeyv(locus.counts.long, cols = the.keys)
        # cbind other samples in, maybe this is very inefficient? maybe cbind whole lot in and exclude the bad rows?
      }
      # I don't think we will eliminate all categories here, so are safe (maybe)
      # Do the chi-sq test
      m2 <- glm(reads ~ bin + affected, data = locus.counts.long[!is.na(affected)], family = "quasipoisson")
      m3 <- glm(reads ~ bin + affected + bin * affected, data = locus.counts.long[!is.na(affected)], family = "quasipoisson")
      test <- anova(m2, m3, test="Chisq")
      
      disease.anova.chisq.statistics[locus.name, sample.name] <- as.numeric(test$Deviance[2])
      disease.anova.chisq.df[locus.name, sample.name]  <- test$Df[2]
      disease.anova.chisq.p.value[locus.name, sample.name]  <- test$Pr[2]
      disease.read.total[locus.name, sample.name] <- locus.counts.long[affected == TRUE, sum(reads)]
    }
    
  }
  str_loglin_test_new(
    data = data,
    loci = loci,
    statistic = disease.anova.chisq.statistics,
    df = disease.anova.chisq.df,
    p.value = disease.anova.chisq.p.value,
    reads.total = disease.read.total, 
    group_case = group_case,
    group_control = group_control,
    group_null = group_null
  )
}

plot.str_loglin_test <- function(x, multi = FALSE, auto.layout = FALSE, 
                                     single.plot = TRUE, # TODO: make this work
                                     statistics = x$p.value, 
                                     diseases = disease.order.by.coverage(x$reads.total), 
                                     read.counts = x$reads.total, 
                                     width = 15, height = 9, 
                                     mfrow = NULL,
                                     mar = c(2.5, 2, 1.5, 1) + 0.1,
                                     plot.blanks = TRUE, 
                                     file = NA, 
                                     read.count.x.weights = c(1, 4),
                                     read.count.y.weights = c(10, 1),
                                     main = paste("Q-Q all loci"),
                                     ...
) {
  assert("x is not of class str_chisq_perm_test", inherits(x, "str_loglin_test"))
  assert("statistics input is NULL. This may be due to $p.value not being defined in data input?", !is.null(statistics))
  assert("diseases input is NULL. This may be due to $reads.total not being defined in data input?", !is.null(diseases))
  if(!is.na(file)) {
    pdf(file, width = width, height = height)
  }
  low.p <- 1 / x$B
  
  if(multi == FALSE) {
    # get all the values together! yo!
    y <- statistics[, x$data$samples[group == x$group_control, sample]]
    if(!is.na(x$group_case)) {
      yy <- statistics[, x$data$samples[group == x$group_case, sample]]
    } else {
      yy <- c()
    }
    qqplot.pvalue(as.vector(as.matrix(y)), pvalues.alt = as.vector(as.matrix(yy)), main = main, ...)
    return()
  }
  pre.par <- par("mfrow")
  if(auto.layout == TRUE && is.null(mfrow)) {
    success <- F
    for(i in 1:10) {
      if(length(diseases) <= i * (2*i + 1)) {
        success <- T
        break
      }
    }
    assert(paste("Too many diseases,", length(diseases), "diseases when the max is 210."), success)
    par(mfrow = c(i, 2 * i + 1))
  }
  if(!is.null(mar)) {
    par(mar = mar)
  }
  max.plots <- prod(par("mfrow"))
  
  plot.count <- 0
  for(disease.name in diseases) {
    
    y <- as.vector(as.matrix(statistics[disease.name, x$data$samples[group == x$group_control, sample]]))
    if(!is.na(x$group_case)) {
      yy <- as.vector(as.matrix(statistics[disease.name, x$data$samples[group == x$group_case, sample]]))
    } else {
      yy <- c()
    }
    
    ## Q-Q plot for Chi^2 data against true theoretical distribution:
    if((plot.count <- plot.count + 1) > max.plots) {
      stop("Error, too many plots attempted")
    }
    qqplot.pvalue(y, pvalues.alt = yy, main = paste(disease.name, "Q-Q"), plot.blanks = plot.blanks)
    if(!is.null(read.counts)) {
      text(weighted.mean(par("usr")[1:2], w = read.count.x.weights), weighted.mean(par("usr")[3:4], w = read.count.y.weights), paste(sum(read.counts[disease.name, ]), "reads", sep="\n"), ...)
    }
  }
  par(mfrow = pre.par)
  if(!is.na(file)) {
    dev.off()
  }
}

