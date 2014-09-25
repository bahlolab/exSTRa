# Functions for analysing repeat expansion data

# Analysing features
# 
require(graphics)
library(ggplot2)
library(xlsx)
library(testit)
library(car)

## functions ##
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

disease.order.by.coverage <- function(X) {
  # X is a data frame of disease coverages, each row is disease
  rownames(X)[order(apply(X, 1, sum), decreasing = T)]
}


plot.different.test.p.values <- function(data.1, data.2, file = NA, width = 15, height = 9) {
  # Plotting the outputs from two different tests together with 
  # computing the Pearson correlation
  if(!is.na(file)) {
    pdf(file, width = width, height = height)
  }
  values.1 <- data.1$p.value
  values.2 <- data.2$p.value
  pre.par <- par("mfrow", "mar")
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
  par(pre.par)
  if(!is.na(file)) {
    dev.off()
  }
}

