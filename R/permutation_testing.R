# Functions to perform the chisq.test permutation testing on strdata objects

str_chisq_permutation_test <- function(data,
                                       diseases = strloci(data)) {
  assert("", inherits(data, "strdata"))
  
  # Pull the input data apart
  diseases <- strdata
  
  #disease.chisq.perm.up.keeprows <- str_chisq_permutation_test_testing(cols = cols.up, keep.rows = keeprows.up, diseases = diseases, features = features, states = states, require.nozero = FALSE)
  disease.chisq.statistics <- data.table()
  disease.chisq.df <- data.table()
  disease.chisq.p.value <- data.table()
  disease.read.total <- data.table()
  for(disease.name in diseases) {
    # Need to first sum the data over the normal samples and the 
    # expanded samples
    #if(disease.name == "FRDA") {
    #  next  
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


str_chisq_permutation_test_testing <- function(cols, keep.rows, diseases, features, states, collapsing.guide = NA, require.nozero = FALSE, B = 100000) {
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
    #  next	
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
