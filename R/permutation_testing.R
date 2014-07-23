# Functions to perform the chisq.test permutation testing on strdata objects

is.str_chisq_perm_test <- function(x) inherits(x, "str_chisq_perm_test")

str_chisq_perm_test_new <- function(data, loci, statistic, df, p.value, reads.total)
{
  assert("Input should be data.frames", is.data.frame(statistic), is.data.frame(df), is.data.frame(p.value), is.data.frame(reads.total))
  structure(
    list(data = data, loci = loci, statistic = statistic, df = df, p.value = p.value, reads.total = reads.total), 
    class = c("str_chisq_perm_test")
  )
}

strloci.str_chisq_perm_test <- function(X) {
  strloci(X$data)
}


str_chisq_permutation_test <- function(data,
                                       cols, 
                                       keep.rows,
                                       loci = strloci(data),
                                       B = 100000, 
                                       require.nozero = TRUE,
                                       group_null = NULL,
                                       group_control = "control",
                                       group_case = "case"
                                       ) {
  # Runs a permutation test
  # Old giving.statistics function inputs:
  # function(cols, keep.rows, diseases, features, states, collapsing.guide = NA, require.nozero = FALSE, B = 100000)

  # Old inputs:
  #disease.chisq.perm.up.keeprows <- str_chisq_permutation_test_testing(cols = cols.up, keep.rows = keeprows.up, diseases = diseases, features = features, states = states, require.nozero = FALSE)

  assert("data is required to be of the class strdata", inherits(data, "strdata"))
  assert("loci should be an atomic vector", is.atomic(loci))
  assert("B should be a single numeric", is.numeric(B), is.atomic(B), length(B) == 1)
  
  # Pull the input data apart
  features <- data$data
  groups <- data$data$group
  # replace states by groups
  
  # Prepare inputs
  loci <- as.character(loci)
  
  # Storing results
  disease.chisq.statistics <- data.table()
  disease.chisq.df <- data.table()
  disease.chisq.p.value <- data.table()
  disease.read.total <- data.table()
  for(locus in loci) {
    # Need to first sum the data over the normal samples and the 
    # expanded samples
    
    # Contingency tables
    cat("STR: ", locus, "\n")
    feature.count <- finding.property.of.features(features, groups, locus, cols, sum) 
    feature.count$sums <- feature.count$result
    
    cont.table <- matrix(feature.count[order(feature.count$group), ]$sums, ncol = nlevels(groups))
    colnames(cont.table) <- unique(sort(feature.count$group))
    rownames(cont.table) <- unique(feature.count$feat)
    
    for(sample.name in levels(features$sample)) {
      #TODO: delete me #sample.data <- subset(features, sample == sample.name & loci == locus) 
      sample.data <- features[.(locus, sample.name)]
      assert("Sample appears to be in multiple groups", length(unique(sample.data$group)) == 1)
      if(is.null(group_null)) {
        # No samples given for null distribution, so use the leave-one-out procedure for the null distribution
        if(sample.data$group[1] == group_case) {
          # expanded sample, compare to normals
          cont.table.1.sample.v <- cont.table[, group_control]
        } else if (sample.data$group[1] == group_control) {
          # normal sample, compare to other normals with leave one out
          leave.one.out.counts <- finding.property.of.features(subset(features, sample != sample.name), as.factor(group_control), locus, cols, sum)
          cont.table.1.sample.v <- leave.one.out.counts[order(leave.one.out.counts$group), ]$result

        } else {
          stop("Unknown group ", sample.data$group[1])
        }
      } else {
        # Have samples for null distribution, so compare controls to the null and cases to the null
        stop("This feature not yet implemented for nulls, controls and cases")
      }
            
      cont.table.1.sample <- matrix(c(unlist(sample.data[, rownames(cont.table), with = F]), cont.table.1.sample.v), ncol = 2)
      colnames(cont.table.1.sample) <- c("subject", "null")
      rownames(cont.table.1.sample) <- rownames(cont.table)

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
      
      disease.chisq.statistics[locus, sample.name] <- test$statistic 
      disease.chisq.df[locus, sample.name]  <- test$parameter
      disease.chisq.p.value[locus, sample.name]  <- test$p.value
      disease.read.total[locus, sample.name] <- sum(cont.table.1.sample[, "subject"])
    }
    
  }
  str_chisq_perm_test_new(
    data = data,
    loci = loci,
    statistic = disease.chisq.statistics,
    df = disease.chisq.df,
    p.value = disease.chisq.p.value,
    reads.total = disease.read.total
  )
}

finding.property.of.features <- function(features, groups, disease.name, cols, FUN) {
  feature.count.fun <- data.frame()
  for (g in levels(groups)) {
    one.sample.one.disease <- subset(features, disease == disease.name & group == g, select = cols)
    column.fun <- c(apply(one.sample.one.disease, 2, FUN))
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