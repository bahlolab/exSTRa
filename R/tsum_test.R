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
#' @import magrittr
#' @export
tsum_test <- function(strscore, 
  trim = 0.15,
  min.quant = 0.5,
  give.pvalue = TRUE, 
  B = 9999, 
  parallel = FALSE, # TRUE for cluster
  cluster = NULL, # a PSOCK cluster created by the parallel package. If NULL, then the
                  # cluster is automatically created.
  cluster_n = NULL # cluster size if cluster == NULL, when NULL, #threads - 1 (but always at least 1)
  ) 
{
  # Check inputs
  assert("strscore should be from class exstra_score.", is.exstra_score(strscore))
  assert("trim should be a number that is at least 0 and less than 0.5.", is.numeric(trim),
    trim >= 0, trim < 0.5)
  assert("min.quant should be a number from 0 to less than 1.", is.numeric(min.quant), min.quant >= 0, min.quant < 1)
  assert("give.pvalue should be logical", is.logical(give.pvalue), !is.na(give.pvalue))
  assert("B should be at least 1 and a whole-number.", B >= 1)
  assert("parallel should be logical.", is.logical(parallel), !is.na(parallel))
  assert("When specified, cluster should be a cluster object from the parallel package.", 
    is.null(cluster) || inherits(cluster, "cluster"))
  assert("cluster_n should be at least 1 and a whole-number.", cluster_n >= 1)
  
  if(parallel) {
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    if(is.null(cluster_n)) {
      # Set the number of cores, max threads - 1 (but at least 1)
      cluster_n <- max(1, n_cores - 1)
    } else {
      if(cluster_n > n_cores) {
        warn.message <- paste0("More threads have been requested by cluster_n (", cluster_n, 
          ") than appears to be available (", 
          n_cores, ").")
        message(warn.message)
      }
    }
  }
  
  # Generate T sum statistic
  T_stats_list <- list()
  for(loc in loci(strscore)) {
    # message("Generating T sum statistics for ", loc)
    qm <- make_quantiles_matrix(strscore, loc = loc, sample = NULL, read_count_quant = 1, 
      method = "quantile7", min.n = 3)
    T_stats_loc <- quant_statistic_sampp(qm, quant = min.quant, trim = trim) # quant at default of 0.5
    T_stats_list[[loc]] <- data.table(sample = names(T_stats_loc), T = T_stats_loc)
  }
  T_stats <- rbindlist(T_stats_list, idcol = "locus")
  
  # Generate p-values (if applicable)
  if(give.pvalue) {
    if(parallel && is.null(cluster)) {
      # create the cluster, just once
      cluster <- makePSOCKcluster(cluster_n)
      cluster_stop <- TRUE
    } else {
      cluster_stop <- FALSE
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
    
    if(parallel && cluster_stop) {
      # stop the cluster that we created
      stopCluster(cluster)
    }
  } else {
    sim.results <- NULL
    pvals <- NULL
  }
  # Prepare output
  exstra_tsum_new_(strscore, T = T_stats, p.values = pvals, 
    qmats = sim.results$qmmats, xecs = sim.results$xecs,
    args = list(trim = trim, min.quant = min.quant, B = B))
}

# test code:
# exp_test <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], B = 50)
#
# exp_test_parallel <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], parallel = TRUE)
#
# exp_test_parallel <- tsum_test(exstra_wgs_pcr_2, parallel = TRUE, B = 50)

# exp_test_all <- tsum_test(exstra_wgs_pcr_2)

# exp_test_parallel_big <- tsum_test(exstra_wgs_pcr_2, parallel = TRUE, B = 9999, cluster_n = 4)
