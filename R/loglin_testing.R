# The log linear tests

is.str_loglin_test <- function(x) inherits(x, "str_loglin_test")

#TODO: make this the actual new function we want
str_loglin_test_new <- function(data, loci, statistic, df, p.value, reads.total,  
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
                                       allow_uncounted_loci = FALSE
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
      sample.data <- locus.counts[sample.name, nomatch = 0] # TODO: this is our problemo
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
      m2 <- glm(reads ~ bin + id, data = locus.counts.long[!is.na(affected)], family = "quasipoisson")
      m3 <- glm(reads ~ bin + id + bin * id, data = locus.counts.long[!is.na(affected)], family = "quasipoisson")
      res <- anova(m2, m3, test="Chisq")
      
      disease.chisq.statistics[locus.name, sample.name] <- test$statistic 
      disease.chisq.df[locus.name, sample.name]  <- test$parameter # TODO: Probably entirely useless??? Maybe replace with remaining categories?
      disease.chisq.p.value[locus.name, sample.name]  <- test$p.value
      disease.read.total[locus.name, sample.name] <- sum(cont.table.1.sample[, "subject"])
    }
    
  }
  str_loglin_test_new(
    data = data,
    loci = loci,
    statistic = disease.chisq.statistics,
    df = disease.chisq.df,
    p.value = disease.chisq.p.value,
    reads.total = disease.read.total, 
    group_case = group_case,
    group_control = group_control,
    group_null = group_null
  )
}