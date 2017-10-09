#' Generate T sum statistics and p-values from simulation. 
#' 
#' When applied to an exstra_score object, T sum statistics are calculated as described
#' in Tankard et al. 
#' May also be applied on a pre-existing exstra_tsum that will regenerate the values.
#' 
#' @import data.table
#' @import stringr
#' @import testit
#' @import parallel
#' @export
tsum_test <- function(strscore, give.pvalue = TRUE, B = 1000, 
  parallel = FALSE, # TRUE for cluster
  cluster = NULL,
  trim = 0.15,
  min.quant = 0.5) 
{
  # Check inputs
  assert("strscore should be from class exstra_score", is.exstra_score(strscore))
  
  
  
  # Generate T sum statistic
  T_stats_list <- list()
  for(loc in loci(strscore)) {
    qm <- make_quantiles_matrix(strscore, loc = loc, sample = NULL, read_count_quant = 1, 
      method = "quantile7", min.n = 3)
    T_stats_loc <- quant_statistic_sampp(qm, quant = quant, trim = trim) # quant at default of 0.5
    T_stats_list[[loc]] <- data.table(sample = names(T_stats_loc), T = T_stats_loc)
  }
  T_stats <- rbindlist(T_stats_list, idcol = "locus")
  
  # Generate p-values (if applicable)
  if(give.pvalue) {
    pvals <- tryCatch(
      simulation_run(multi_cohorts_combinations[[this.cohort]], 
        B = B, 
        subject_in_background = TRUE, 
        sort_sim_qm = TRUE, 
        use_truncated_sd = FALSE,
        truncated_sd_fudge = NA,
        use_truncated_sd_in_quant_statistic = FALSE,
        median_mad = TRUE, 
        parallel = TRUE, # TRUE for cluster
        cluster = cl,
        trim = trim
      ),
      error = function(x) { list(cohort = this.cohort, p.matrix = c(0, 1), error = x) }
    )
  } else {
    pvals <- NULL
  }
  # Prepare output
  exstra_tsum_new_(strscore, T = T_stats, p = pvals)
}

# test code:
# tsum_test(exstra_wgs_pcr_2["HD"], give.pvalue = FALSE)


