# Functions for analysing repeat expansion data

# Analysing features
# 
require(graphics)
library(ggplot2)
library(xlsx)
library(testit)
library(car)

## functions ##
finding.property.of.features <- function(features, states, disease.name, cols, FUN) {
  feature.count.fun <- data.frame()
  for (stay in states) {
    one.sample.one.disease <- subset(features, disease == disease.name & state == stay, select = cols)
    column.fun <- c(apply(one.sample.one.disease, 2, FUN))
    feature.count.fun <- rbind(feature.count.fun, column.fun) 
  }
  names( feature.count.fun ) <- cols
  
  feature.count.fun.for.plot <- data.frame()
  for(i in 1:dim(feature.count.fun)[2]) {
    a <- data.frame( feat = rep(names(feature.count.fun)[i], 2), result = feature.count.fun[,i], state = states)
    feature.count.fun.for.plot <- rbind(feature.count.fun.for.plot, a)
  }
  return(feature.count.fun.for.plot)
}

collapse.contigency.table.2 <- function(c.table, rows, names = c("provided", "other")) {
  # collapse out a table so there are only two rows
  # The given row/s are one group, with all the others the other group
  if(! is.numeric(rows)) {
    rows <- which(is.element(rownames(c.table), rows))
  }
  out <- rbind(apply(c.table[rows, , drop=FALSE], 2, sum), apply(c.table[-rows, , drop=FALSE], 2, sum))
  rownames(out) <- names
  out
}

collapse.contigency.table.single <- function(c.table, rows, name = "") {
  # collapse out a table so that the given rows are collapsed, leaving others in tact
  # The given row/s are one group, with all the others the other group
  if(! is.numeric(rows)) {
    rows <- which(is.element(rownames(c.table), rows))
  }
  out <- rbind(apply(c.table[rows, , drop=FALSE], 2, sum), c.table[-rows, , drop=FALSE])
  rownames(out)[1] <- name
  out
}

remove.contigency.table.rows <- function(c.table, rows) {
  # remove the given rows, leaving others in tact
  # The given row/s are one group, with all the others the other group
  if(! is.numeric(rows)) {
    rows <- which(is.element(rownames(c.table), rows))
  }
  c.table[-rows, ]
}

simplify.bases <- function(x) {
  x <- as.character(x)
  if(nchar(x) > 4) {
    paste0(nchar(x), "bp")
  } else {
    x  
  }
}


giving.statistics <- function(cols, keep.rows, diseases, features, states, collapsing.guide = NA, require.nozero = FALSE, B = 100000) {
  pair.type <- NA  
  # Function for permutation testing
  #cat("Performing", pair.type, "analysis.\n")
  set.seed(c(32,3461,334345,982,0209,91859))
  disease.chisq.statistics <- data.frame()
  disease.chisq.df <- data.frame()
  disease.chisq.p.value <- data.frame()
  disease.read.total <- data.frame()
  for(disease.name in diseases) {
    # Need to first sum the data over the normal samples and the 
    # expanded samples
    #if(disease.name == "FRDA") {
    #	next	
    #}
    # Contingency tables
    cat("STR: ", disease.name, "\n")
    feature.count.for.plot <- finding.property.of.features(features, states, disease.name, cols, sum) 
    feature.count.for.plot$sums <- feature.count.for.plot$result
    
    cont.table <- matrix(feature.count.for.plot[order(feature.count.for.plot$state), ]$sums, ncol = 2)
    colnames(cont.table) <- unique(sort(feature.count.for.plot$state))
    rownames(cont.table) <- unique(feature.count.for.plot$feat)
    
    #if(!is.na(collapsing.guide[1])) {
    #	for(simple_feature in names(collapsing.guide)) {
    #		cont.table <- collapse.contigency.table.single(cont.table, 
    #			grep(collapsing.guide[simple_feature], rownames(cont.table), value = T), name = simple_feature)
    #	}
    #}
    # TODO:
    # Loop over the 20 normal and 20 affected samples and see if we can get significant statistics
    if(!is.null(features$disease)) {
      str.names <- features$disease
    } else if(!is.null(features$STR)) {
      str.names <- features$STR
    } else {
      stop("Could not find the str column in the features of the giving.statistics function.")
    }
    for(sample.name in levels(features$sample)) {
      sample.data <- subset(features, sample == sample.name & str.names == disease.name)
      if(sample.data$state[1] == "expanded") {
        # expanded sample, compare to normals
        cont.table.1.sample <- cont.table
      } else {
        # normal sample, compare to other normals with leave one out
        leave.one.out.counts <- finding.property.of.features(subset(features, sample != sample.name), states, disease.name, cols, sum)
        cont.table.1.sample <- matrix(leave.one.out.counts[order(leave.one.out.counts$state), ]$result, ncol = 2)
        colnames(cont.table.1.sample) <- unique(sort(leave.one.out.counts$state))
        rownames(cont.table.1.sample) <- unique(leave.one.out.counts$feat)
      }		
      colnames(cont.table.1.sample)[1] <- "subject"
      cont.table.1.sample[, "subject"] <- unlist(sample.data[rownames(cont.table.1.sample)])
      if(!is.na(collapsing.guide[1])) {
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
      
      disease.chisq.statistics[disease.name, sample.name] <- test$statistic 
      disease.chisq.df[disease.name, sample.name]  <- test$parameter
      disease.chisq.p.value[disease.name, sample.name]  <- test$p.value
      disease.read.total[disease.name, sample.name] <- sum(cont.table.1.sample[, "subject"])
    }
    
  }
  list(statistic = disease.chisq.statistics, df = disease.chisq.df, p.value = disease.chisq.p.value, reads.total = disease.read.total)
}

disease.order.by.coverage <- function(X) {
  # X is a data frame of disease coverages, each row is disease
  rownames(X)[order(apply(X, 1, sum), decreasing = T)]
}

plot.qq.diseases <- function(file = NA, 
                             data = NULL, 
                             statistics = data$p.value, 
                             diseases = disease.order.by.coverage(data$reads.total), 
                             read.counts = data$reads.total, 
                             q.fun, 
                             width = 15, height = 9, 
                             mfrow=c(3, 7), 
                             mar = c(2.5, 2, 1.5, 1) + 0.1, 
                             log = "", 
                             plot.blanks = TRUE, 
                             xlim = NULL, 
                             ylim = NULL,
                             low.p = 1 / 100000
) {
  assert("statistics input is NULL. This may be due to $p.value not being defined in data input?", !is.null(statistics))
  assert("diseases input is NULL. This may be due to $reads.total not being defined in data input?", !is.null(diseases))
  if(!is.na(file)) {
    pdf(file, width = width, height = height)
  }
  par(mfrow = mfrow, mar = mar)
  plot.count <- 0
  for(disease.name in diseases) {
    #if(disease.name == "FRDA") {
    #	plot(NA, NA, main = "FRDA")
    #	next	
    #}
    #cat("working on", disease.name, "\n")
    if(sum(grepl("expanded|normal", colnames(statistics))) == dim(statistics)[2]) {
      d <- data.frame(
        expanded = unlist(statistics[disease.name, grepl("expanded", colnames(statistics))]), 
        normal = unlist(statistics[disease.name, grepl("normal", colnames(statistics))])
      )
      y <- d$normal
      yy <- d$expanded
    } else {
      d <- data.frame(
        samples = unlist(statistics[disease.name, colnames(statistics)])
      )
      y <- d$samples
      yy <- numeric()
    }
    
    ## Q-Q plot for Chi^2 data against true theoretical distribution:
    axis.limits <- c(0, -log10(low.p))
    if(is.null(xlim)) { xlim <- axis.limits }
    if(is.null(ylim)) { ylim <- axis.limits }
    if(sum(is.na(y)) >= length(y)) {
      if(plot.blanks) {
        plot(NA, xlim = xlim, ylim = ylim, main = paste(disease.name, "Q-Q"), log = log)
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black")
      }
    } else {
      if((plot.count <- plot.count + 1) > prod(mfrow)) {
        stop("Error, too many plots attempted")
      }
      #qqplot(-log10(q.fun(ppoints(length(y)))), 
      #       -log10(y),
      #       main = paste(disease.name, "Q-Q")
      #       , xlim = xlim
      #       , ylim = ylim
      #)
      #curve(x * 1, 0, -log10(low.p), add = TRUE, col = "red")
      qqplot.pvalue(y, pvalues.alt = yy, main = paste(disease.name, "Q-Q"))
      if(!is.null(read.counts)) {
        text(.154, 4, paste(sum(read.counts[disease.name, ]), "reads", sep="\n"), cex = 2)
      }
    }
    #if(sum(is.na(yy)) < length(yy)) {
    #  points(-log10(q.fun(ppoints(length(yy)))), -log10(sort(yy)), col = "blue", pch = 2)
    #}
  }
  if(!is.na(file)) {
    dev.off()
  }
}

plot.differet.test.p.values <- function(data.1, data.2, file = NA, width = 15, height = 9) {
  # Plotting the outputs from two different tests together with 
  # computing the Pearson correlation
  if(!is.na(file)) {
    pdf(file, width = width, height = height)
  }
  values.1 <- data.1$p.value
  values.2 <- data.2$p.value
  par(mfrow=c(3, 7), mar = c(5, 4, 4, 2) + 0.1)
  lims <- function(X) {
    if(sum(!is.na(X)) == 0) {
      c(1e-5, 1)
    } else {
      range(X, na.rm = T)  
    }
  }
  for(disease in diseases.normal.size.ordered) {
    x.1 <- unlist(values.1[disease, ])
    x.2 <- unlist(values.2[disease, ])
    correla <- cor(x.1, x.2, use = "na.or.complete")
    plot(x.1, x.2, log = "xy", main = paste0(disease, "\nR=", sprintf("%0.4f", correla)), xlab = "", ylab = "", xlim = lims(x.1), ylim = lims(x.2))
  }
  if(!is.na(file)) {
    dev.off()
  }
}


qqplot.pvalue <- function(pvalues = NULL, p.trans = -log10(pvalues), 
                          pvalues.alt = NA, p.trans.alt = -log10(pvalues.alt), 
                          main = "QQ plot of p-values", xlab = "Null", ylab = "Observed") {
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
  plot(NULL, ylim=c(0, MAX.y), xlim=c(0, MAX.x), axes=FALSE, xlab="", ylab="")
  grid()
  lines(null, -log10(c95), col = "red", lty = 2)
  lines(null, -log10(c05), col = "red", lty = 2)
  
  ## add the diagonal
  abline(0,1,col="red")
  par(new=T)
  
  ## add the qqplot
  qqplot(null, p.trans, ylim=c(0, MAX.y), xlim=c(0, MAX.x), main=main, xlab = xlab, ylab = ylab)
  
  # Add the alternative QQ plot if given
  if(sum(is.na(p.trans.alt)) < length(p.trans.alt)) {
    points(-log10(qunif(ppoints(length(p.trans.alt)))), sort(p.trans.alt, decreasing = TRUE), col = "blue", pch = 2)
  }
}
