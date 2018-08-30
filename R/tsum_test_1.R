#' Generate T sum statistics and p-values from simulation. 
#' 
#' Original thesis version. 
#' When applied to an exstra_score object, T sum statistics are calculated as described
#' in Tankard et al. 
#' May also be applied on a pre-existing exstra_tsum that will regenerate the values.
#' 
#' @param strscore An exstra_score object.
#' @param trim Trim this proportion of data points at each quantile 
#'             level (rounded up). Must be at least 0 and less than 0.5, but values close
#'             to 0.5 may remove all samples and hence result in an error.  
#' @param trim.cc Trim value used in case-control analysis. Default of 0.
#' @param min.quant Only quantiles above this value are used in constructing the statistic.  
#' @param give.pvalue Whether to calculate the p-value. As this can be slow it can be
#'                    useful to turn off if only the t sum statistics are required. 
#' @param B Number of simulations in calculating null distributions. The denominator will
#'          be B + 1, hence values of B = 10^i - 1 will result in p-values that are 
#'          decimal fractions. 
#' @param correction Correction method of p_value() function.
#' @param alpha Signficance level of p_value() function.
#' @param case_control If TRUE, only calculate for samples designated cases. Otherwise
#'                 all samples are used to calculate the background distribution.
#' @param parallel Use the parallel package when simulating the distribution, creating the
#'                 required cluster. 
#'                 If cluster is specified then this option makes no difference. 
#' @param cluster_n If parallel is TRUE, then the number of nodes in the cluster is 
#'                  automatically set as 1 less than those available on your machine. 
#'                  (but never less than 1). This option allows manual setting of the 
#'                  number of nodes, either less to free up other resources, or more to 
#'                  maximize available resources. 
#'                 If cluster is specified then this option makes no difference.
#' @param cluster  A cluster object from the parallel package. Use if you wish to set up 
#'                 the cluster yourself or reuse an existing cluster. 
#' @param keep.sim.tsum For inspection of simulations. 
#'                      If TRUE, keep all simulation Tsum statistics in output$xecs (default FALSE).
#'                      
#' @return An exstra_tsum object with T statistics and p-values (if calculated).
#' 
#' @seealso \code{\link{tsum_p_value_summary}}
#' 
#' @references 
#'         Rick M. Tankard, Martin B. Delatycki, Paul J. Lockhart, 
#'         Melanie Bahlo. 
#'         Detecting known repeat expansions with standard protocol next generation 
#'         sequencing, towards developing a single screening test for neurological repeat 
#'         expansion disorders. 
#'         \emph{bioRxiv} 157792; 
#'         doi: \url{https://doi.org/10.1101/157792}
#' 
#' @examples 
#' exp_test <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], B = 50)
#' exp_test
#'
#' \dontrun{
#' exp_test_parallel <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], parallel = TRUE, B = 999)
#' exp_test_parallel
#' }
#' 
#' @import data.table
#' @import stringr
#' @import testit
#' @import parallel
#' 
tsum_test_1 <- function(strscore, 
  trim = 0.15,
  trim.cc = 0,
  min.quant = 0.5,
  give.pvalue = TRUE, 
  B = 999, 
  correction = c("bf", "loci", "samples", "uncorrected"),
  alpha = 0.05,
  case_control = FALSE, 
  parallel = FALSE, # TRUE for cluster
  cluster_n = NULL, # Cluster size if cluster == NULL. When NULL, #threads - 1 (but always at least 1)
  cluster = NULL, # As created by the parallel package. If cluster == NULL and parallel == TRUE, then a
                  # PSOCK cluster is automatically created with the parallel package.
  keep.sim.tsum = FALSE
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
  assert("cluster_n should be at least 1 and a whole-number.", is.null(cluster_n) || cluster_n >= 1)
  
  # temp code until this is developed
  if(case_control && give.pvalue) {
    stop("tsum_test() cannot yet use cases and controls to generate p-values from simulation.") 
  }
  
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
  
  if(inherits(cluster, "cluster")) {
    # as a cluster has been given, assume we actually do want to use the parallel package
    parallel <- TRUE
  }
  
  # Generate T sum statistic
  T_stats_list <- list()
  for(loc in loci(strscore)) {
    # message("Generating T sum statistics for ", loc)
    qm <- make_quantiles_matrix(strscore, loc = loc, sample = NULL, read_count_quant = 1, 
      method = "quantile7", min.n = 3)
    if(case_control) {
      # Only calculating for case samples
      T_stats_loc <- quant_statistic_sampp(qm, quant = min.quant, trim = trim.cc,
        case_samples = strscore$samples[group == "case", sample]) 
    } else {
      # Using all samples as the background:
      T_stats_loc <- quant_statistic_sampp(qm, quant = min.quant, trim = trim) # quant at default of 0.5
    }
    T_stats_list[[loc]] <- data.table(sample = names(T_stats_loc), tsum = T_stats_loc)
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
        trim = trim,
        T_stats = T_stats
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
  outtsum <- exstra_tsum_new_(strscore, tsum = T_stats, p.values = pvals, 
    qmats = sim.results$qmmats, xecs = sim.results$xecs,
    correction = correction,
    alpha = alpha, 
    args = list(trim = trim, min.quant = min.quant, B = B))
  if(give.pvalue) {
    Nsim <- B * strscore$samples[, .N]
    outtsum$stats[, p.value.sd := 
        sqrt(p.value * ((Nsim + 2)/(Nsim + 1) - p.value) / Nsim)
      ]
  }
  if(! keep.sim.tsum) {
    for(i in seq_along(outtsum$xecs)) {
      outtsum$xecs[[i]]$sim_T <- NULL
    }
  }
  outtsum
}

# test code:
# exp_test <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], B = 50)
#
# exp_test_parallel <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], parallel = TRUE)
#
# exp_test_parallel <- tsum_test(exstra_wgs_pcr_2, parallel = TRUE, B = 50)

# exp_test_all <- tsum_test(exstra_wgs_pcr_2)

# exp_test_parallel_big <- tsum_test(exstra_wgs_pcr_2, parallel = TRUE, B = 9999, cluster_n = 4)
