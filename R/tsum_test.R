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
tsum_test <- function(strscore, 
  trim = 0.15,
  min.quant = 0.5,
  give.pvalue = TRUE, 
  B = 1000, 
  parallel = FALSE, # TRUE for cluster
  cluster = NULL,
  cluster_n = NULL # cluster size if cluster == NULL, when NULL, #threads - 1 (but always at least 1)
  ) 
{
  # Check inputs
  assert("strscore should be from class exstra_score", is.exstra_score(strscore))
  
  
  
  # Generate T sum statistic
  T_stats_list <- list()
  for(loc in loci(strscore)) {
    qm <- make_quantiles_matrix(strscore, loc = loc, sample = NULL, read_count_quant = 1, 
      method = "quantile7", min.n = 3)
    T_stats_loc <- quant_statistic_sampp(qm, quant = min.quant, trim = trim) # quant at default of 0.5
    T_stats_list[[loc]] <- data.table(sample = names(T_stats_loc), T = T_stats_loc)
  }
  T_stats <- rbindlist(T_stats_list, idcol = "locus")
  
  # Generate p-values (if applicable)
  if(give.pvalue) {
    if(parallel && is.null(cluster_n)) {
      cluster_n <- max(1, detectCores(all.tests = FALSE, logical = TRUE) - 1)
    }
    sim.results <- tryCatch(
      simulation_run(strscore, 
        B = B, 
        subject_in_background = TRUE, 
        sort_sim_qm = TRUE, 
        use_truncated_sd = FALSE,
        truncated_sd_fudge = NA,
        use_truncated_sd_in_quant_statistic = FALSE,
        median_mad = TRUE, 
        parallel = parallel, # TRUE for cluster
        cluster = cluster,
        cluster_n = cluster_n,
        trim = trim
      ),
      error = function(x) { 
          # list(cohort = this.cohort, p.matrix = c(0, 1), error = x)
          warning(x)
          return(NULL)
        }
    )
    pvals <- sim.results$p.matrix
  } else {
    sim.results <- NULL
    pvals <- NULL
  }
  # Prepare output
  exstra_tsum_new_(strscore, T = T_stats, p = pvals, 
    qmats = sim.results$qmmats, xecs = sim.results$xecs)
}

# test code:
# tsum_test(exstra_wgs_pcr_2["HD"], give.pvalue = FALSE)


