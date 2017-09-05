# Functions likely of no use, poor performance


#
exstra_score_ks_tests <- function(rsc, locus = NULL, controls = c("control", "all")) {
  # Performs Kolmogorov-Smirnov Tests on samples, comparing to other samples
  #
  # controls allows either just control samples to be used as the population distribution,
  # or all other samples including designated cases. This makes no difference if there is
  # only a single case. 
  if(!is.exstra_score(rsc)) {
    stop("rsc is not object of class exstra_score")
  }
  if(is.null(locus)) {
    strlocis <- loci(rsc)
  } else {
    strlocis <- locus
  } 
  results <- data.table(
    locus = rep(strlocis, length(rsc$samples[group == 'case', sample])), 
    sample = rep(rsc$samples[group == 'case', sample], each = length(strlocis)), 
    p.value = NA_real_ #,
    #test = list(list())
  )
  setkey(results, locus, sample)
  for(loc in strlocis) {
    loc_scores <- rsc$data[locus == loc]
    for(samp in rsc$samples[group == 'case', sample]) {
      KS <- ks.test(loc_scores[group == "control"]$rep, loc_scores[sample == samp]$rep, 
        alternative = "greater", exact = NULL)
      results[list(loc, samp), p.value := KS$p.value]
      #results[list(loc, samp), test := KS]
    }
  }
  
  return(
    data.frame(results)
  )
}