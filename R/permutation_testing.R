# Functions to perform the chisq.test permutation testing on strdata objects

#TODO: use locus instead of disease for location name (maybe change all instances of disease to locus?)

is.str_chisq_perm_test <- function(x) inherits(x, "str_chisq_perm_test")

str_chisq_perm_test_new <- function(data, loci, statistic, df, p.value, reads.total, B, 
                                    group_case,
                                    group_control,
                                    group_null)
{
  assert("Input should be data.frames", is.data.frame(statistic), is.data.frame(df), is.data.frame(p.value), is.data.frame(reads.total))
  assert("B should be a single positive numeric", is.numeric(B), is.atomic(B), length(B) == 1, B >= 1)
  assert("data should be of class strdata", is.strdata(data))
  assert("loci should be a character vector", is.character(loci))
  structure(
    list(data = data, loci = loci, statistic = statistic, df = df, p.value = p.value, reads.total = reads.total, B = B, 
         group_case = group_case,
         group_control = group_control,
         group_null = group_null), 
    class = c("str_chisq_perm_test")
  )
}

strloci.str_chisq_perm_test <- function(X) {
  strloci(X$data)
}


str_chisq_permutation_test <- function(data,
                                       cols, 
                                       keep.rows,
                                       collapsing.guide = NULL, 
                                       loci = strloci(data), 
                                       B = 100000, 
                                       require.nozero = TRUE,
                                       group_null = NULL,
                                       group_control = "control",
                                       group_case = "case",
                                       allow_uncounted_loci = FALSE
                                       ) {
  # Runs a permutation test
  # Old giving.statistics function inputs:
  # function(cols, keep.rows, diseases, features, states, collapsing.guide = NA, require.nozero = FALSE, B = 100000)

  # Old inputs:
  #disease.chisq.perm.up.keeprows <- str_chisq_permutation_test_testing(cols = cols.up, keep.rows = keeprows.up, diseases = diseases, features = features, states = states, require.nozero = FALSE)

  assert("data is required to be of the class strdata", inherits(data, "strdata"))
  assert("loci should be an atomic vector", is.atomic(loci))
  assert("B should be a single positive numeric", is.numeric(B), is.atomic(B), length(B) == 1, B >= 1)

  
  # Pull the input data apart
  features <- data$data
  groups <- data$data$group
  # replace states by groups
  
  # Prepare inputs
  loci <- as.character(loci)
  B <- as.integer(B)
  
  # Warnings
  if(length(setdiff(loci, levels(features$locus))) > 0) {
    mess <- "loci are given that are not in the data"
    if(allow_uncounted_loci) {
      warning(mess)      
    } else {
      stop(mess)
    }
  }
  
  # Storing results
  disease.chisq.statistics <- data.frame() #TODO: make data.table
  disease.chisq.df <- data.frame()
  disease.chisq.p.value <- data.frame()
  disease.read.total <- data.frame()
  for(locus.name in loci) {
    # Need to first sum the data over the normal samples and the 
    # expanded samples
    
    # Contingency tables
    cat("STR: ", locus.name, "\n")
    feature.count <- finding.property.of.features(features, groups, locus.name, cols, sum) 
    feature.count$sums <- feature.count$result
    
    cont.table <- matrix(feature.count[order(feature.count$group), ]$sums, ncol = nlevels(groups))
    colnames(cont.table) <- unique(sort(feature.count$group))
    rownames(cont.table) <- unique(feature.count$feat)
    
    for(sample.name in levels(features$sample)) {
      sample.data <- features[.(locus.name, sample.name), nomatch = 0] # TODO: this is our problemo
      if(dim(sample.data)[1] == 0) {
        next
      }
      assert("Sample appears to be in multiple groups", length(unique(sample.data$group)) == 1)
      if(is.null(group_null)) {
        # No samples given for null distribution, so use the leave-one-out procedure for the null distribution
        if(sample.data$group[1] == group_case) {
          # expanded sample, compare to normals
          cont.table.1.sample.v <- cont.table[, group_control]
        } else if (sample.data$group[1] == group_control) {
          # normal sample, compare to other normals with leave one out
          leave.one.out.counts <- finding.property.of.features(subset(features, sample != sample.name), as.factor(group_control), locus.name, cols, sum)
          cont.table.1.sample.v <- leave.one.out.counts[order(leave.one.out.counts$group), ]$result

        } else {
          stop("Unknown group ", sample.data$group[1])
        }
      } else {
        # Have samples for null distribution, so compare controls to the null and cases to the null
        # set: cont.table.1.sample.v
        stop("This feature not yet implemented for nulls, controls and cases")
      }
            
      #cont.table.1.sample <- matrix(c(unlist(sample.data[, rownames(cont.table), with = F]), cont.table.1.sample.v), ncol = 2)
      cont.table.1.sample <- matrix(c(sapply(sample.data[, rownames(cont.table), with = F], sum), cont.table.1.sample.v), ncol = 2)
      
      colnames(cont.table.1.sample) <- c("subject", "null")
      rownames(cont.table.1.sample) <- rownames(cont.table)

      if(!is.null(collapsing.guide)) {
        for(simple_feature in names(collapsing.guide)) {
          cont.table.1.sample <- collapse.contigency.table.single(cont.table.1.sample, 
                                                                  grep(collapsing.guide[simple_feature], rownames(cont.table.1.sample), value = T), name = simple_feature)
        }
      }
      cont.table.1.sample <- cont.table.1.sample[keep.rows, , drop = F]
      if(require.nozero) {
        cont.table.1.sample <- cont.table.1.sample[apply(cont.table.1.sample, 1, sum) != 0, , drop = F] # Remove zero rows
      }
      if(dim(cont.table.1.sample)[1] == 0 || sum(cont.table.1.sample[, 1]) == 0) {
        # No rows left
        test <- 
          list(  statistic = NA
                 , parameter = NA
                 , p.value = NA
          )
      } else {
        # Do the chi-sq test
        test <- chisq.test(cont.table.1.sample, simulate.p.value = T, B = B)
      }
      
      disease.chisq.statistics[locus.name, sample.name] <- test$statistic 
      disease.chisq.df[locus.name, sample.name]  <- test$parameter # TODO: Probably entirely useless??? Maybe replace with remaining categories?
      disease.chisq.p.value[locus.name, sample.name]  <- test$p.value
      disease.read.total[locus.name, sample.name] <- sum(cont.table.1.sample[, "subject"])
    }
    
  }
  str_chisq_perm_test_new(
    data = data,
    loci = loci,
    statistic = disease.chisq.statistics,
    df = disease.chisq.df,
    p.value = disease.chisq.p.value,
    reads.total = disease.read.total, 
    B = B, 
    group_case = group_case,
    group_control = group_control,
    group_null = group_null
  )
}

finding.property.of.features <- function(features, groups, locus.name, cols, FUN) {
  feature.count.fun <- data.frame()
  for (g in levels(groups)) {
    one.sample.one.locus <- subset(features, locus == locus.name & group == g, select = cols)
    column.fun <- c(apply(one.sample.one.locus, 2, FUN))
    feature.count.fun <- rbind(feature.count.fun, column.fun) 
  }
  names( feature.count.fun ) <- cols
  
  feature.count.fun.for.plot <- data.frame()
  for(i in 1:dim(feature.count.fun)[2]) {
    a <- data.frame( feat = rep(names(feature.count.fun)[i], nlevels(groups)), result = feature.count.fun[,i], group = levels(groups))
    feature.count.fun.for.plot <- rbind(feature.count.fun.for.plot, a)
  }
  return(feature.count.fun.for.plot)
}



qqplot.pvalue <- function(pvalues = NULL, p.trans = -log10(pvalues), 
                          pvalues.alt = NA, p.trans.alt = -log10(pvalues.alt), 
                          main = "QQ plot of p-values", xlab = "Null", ylab = "Observed", 
                          plot.blanks = FALSE, ...) {
  # Create a QQ plot of p-values with confidence intervals
  # http://gettinggeneticsdone.blogspot.com.au/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
  # 17th July 2014
  N <- length(p.trans) ## number of p-values
  assert("Require at least one observation", N > 0)
  
  ## create the null distribution 
  ## (-log10 of the uniform)
  null <- -log10(ppoints(N))
  MAX.x <- max(null)
  MAX.y <- max(MAX.x, p.trans, p.trans.alt, na.rm = T)
  
  ## create the confidence intervals
  c95 <- rep(0,N)
  c05 <- rep(0,N)
  
  ## the jth order statistic from a 
  ## uniform(0,1) sample 
  ## has a beta(j,n-j+1) distribution 
  ## (Casella & Berger, 2002, 
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:N){
    c95[i] <- qbeta(0.95,i,N-i+1)
    c05[i] <- qbeta(0.05,i,N-i+1)
  }
  
  ## plot the two confidence lines
  plot(NULL, ylim=c(0, MAX.y), xlim=c(0, MAX.x), axes=FALSE, xlab="", ylab="", ...)
  grid()
  lines(null, -log10(c95), col = "red", lty = 2)
  lines(null, -log10(c05), col = "red", lty = 2)
  
  ## add the diagonal
  abline(0,1,col="red")
  par(new=T)
  
  ## add the qqplot
  qqplot(null, p.trans, ylim=c(0, MAX.y), xlim=c(0, MAX.x), main=main, xlab = xlab, ylab = ylab, ...)
  
  # Add the alternative QQ plot if given
  if(sum(is.na(p.trans.alt)) < length(p.trans.alt)) {
    points(-log10(qunif(ppoints(length(p.trans.alt)))), sort(p.trans.alt, decreasing = TRUE), col = "blue", pch = 2)
  }
}

plot.str_chisq_perm_test <- function(x, multi = FALSE, auto.layout = FALSE, 
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
  assert("x is not of class str_chisq_perm_test", inherits(strcount.perm, "str_chisq_perm_test"))
  assert("statistics input is NULL. This may be due to $p.value not being defined in data input?", !is.null(statistics))
  assert("diseases input is NULL. This may be due to $reads.total not being defined in data input?", !is.null(diseases))
  if(!is.na(file)) {
    pdf(file, width = width, height = height)
  }
  low.p <- 1 / x$B

  if(multi == FALSE) {
    # get all the values together! yo!
    #strcount.perm$data$data[group == "control", c("sample", "disease", "group"), with = F]
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
