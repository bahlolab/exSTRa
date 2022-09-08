#' Generate T sum statistics and p-values from simulation. 
#' 
#' When applied to an exstra_score object, T sum statistics are calculated as described
#' in Tankard et al. 
#' May also be applied on a pre-existing exstra_tsum that will regenerate the values.
#' 
#' @param strscore An exstra_score object.
#' @param trim Trim this proportion of data points at each quantile 
#'             level (rounded up). Must be at least 0 and less than 0.5, but values close
#'             to 0.5 may remove all samples and hence result in an error.  
#' @param trim.cc Trim value used in case-control analysis. Default of 0.
#' @param trim.all Trim value when looking for outliers in all samples. Default 0.15.
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
#' @param early_stop Simulation may use less replicates when all p-values are large, controlled with early_A.
#' @param early_A  Simulations may stop when p.value.sd < early_A * min(p.value). 
#'                 Checked approximately when the number of simulations has doubled.
#' @param early_stop_min Minimum number of simulations to run before early termination.
#' @param parallel Use the parallel package when simulating the distribution, creating the
#'                 required cluster. 
#'                 If cluster is specified then this option makes no difference. 
#' @param cluster_n If parallel is TRUE, then the number of nodes in the cluster is 
#'                  automatically set as half of those available on your machine 
#'                  (but never less than 1). This option allows manual setting of the 
#'                  number of nodes, either less to free up other resources, or to 
#'                  maximize available resources. 
#'                 If cluster is specified then this option makes no difference.
#' @param cluster A snow cluster created by snow::makeCluster(). Use if you wish to set up 
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
#' @import snow
#' @export
tsum_test <- function(strscore, 
                           trim = ifelse(case_control, trim.cc, trim.all),
                           trim.all = 0.15,
                           trim.cc = 0,
                           min.quant = 0.5,
                           give.pvalue = TRUE, 
                           B = 999, 
                           correction = c("bf", "loci", "samples", "uncorrected"),
                           alpha = 0.05,
                           case_control = FALSE,
                           early_stop = TRUE,
                           early_A = 0.25,
                           early_stop_min = 50,
                           parallel = FALSE, # TRUE for cluster
                           cluster_n = NULL, # Cluster size if cluster == NULL. When NULL, #threads / 2 (but always at least 1)
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
  
  # Warning for parallel usage
  if(parallel) {
    warning("Use of tsum_test(parallel = TRUE) may not be beneficial.\n",
            "  Optimizations in serial code have reduced the benefits of parallization.\n",
            "  Additionally, often R jobs tend to be idle which has not been resolved.",
            immediate. = TRUE)
  }
  
  # trim too high?
  if (trim > 0.3) {
    #TODO: make this a better test, such as the number of samples left (maybe I've already done this?) - Rick
    warning("trim is set to ", trim, ", removing ", trim * 200, "% of the data.")
  }
  if(case_control) {
    n_back_samples <- length(trim_vector(strscore$samples[group != "case", .N], trim))
  } else {
    n_back_samples <- length(trim_vector(strscore$samples[, .N], trim))
  }
  if(n_back_samples < 6) {
    warning("Trimming at each quantile only leaves ", n_back_samples, " observations at each level.")
  }
  
  # Check samples after cases are removed in case_control
  if(case_control && strscore$samples[group != "case", .N] < 5) {
    warning("May have too few control samples (n = ", strscore$samples[group != "case", .N], ").\n",
            "Same tsum statistics may be NaN.")
  }
  
  # Set or check the number of cores in parallel, when no cluster is specified
  if(inherits(cluster, "cluster")) {
    # as a cluster has been given, assume we actually do want to use the parallel package
    parallel <- TRUE
  } else {
    if(parallel) {
      n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
      if(is.null(cluster_n)) {
        # Set the number of cores, max threads / 2 (but at least 1)
        cluster_n <- max(1, n_cores / 2)
      } else {
        if(cluster_n > n_cores) {
          warn.message <- paste0("More threads have been requested by cluster_n (", cluster_n, 
                                 ") than appears to be available (", 
                                 n_cores, ").")
          message(warn.message)
        }
      }
    }
  }
  
  if(give.pvalue && parallel) {
    # Create a new PSOCKcluster if required
    if(is.null(cluster)) {
      # create the cluster, just once
      cluster <- snow::makeCluster(cluster_n)
      on.exit(stopCluster(cluster), add = TRUE)
    }
    
    # Load required functions onto cluster nodes 
    # TODO: is this useless, since this doesn't bring functions into namespace?
    clusterEvalQ(cluster, { 
      requireNamespace(testit); 
      requireNamespace(magrittr);
      requireNamespace(data.table);
      requireNamespace(exSTRa) 
    })
  }
  
  # Generate T sum statistic
  if(parallel) {
    all_loci <- loci(strscore)
    names(all_loci) <- all_loci
    # Save copying the whole object each time
    strscore_loc_list <- lapply(all_loci, function(loc) { strscore[loc]} )
    tsum_statistic_1locus_failsafe <- function(...) {
      tryCatch(tsum_statistic_1locus(...), error = function(e) { NULL })
    }
    T_stats_list <- snow::parLapply(cl = cluster,
                                    strscore_loc_list, 
                                    tsum_statistic_1locus_failsafe,
                                    min.quant = min.quant,
                                    case_control = case_control, trim = trim,
                                    give.pvalue = give.pvalue, B = B,
                                    early_stop = early_stop, early_A = early_A, min_stop = early_stop_min,
                                    verbose = FALSE)
    
    # Remove errored runs
    T_stats_list <- T_stats_list[!sapply(T_stats_list, is.null)]
    
  } else {
    T_stats_list <- vector('list', length(loci(strscore)))
    names(T_stats_list) <- loci(strscore)
    for(loc in loci(strscore)) {
      message("Working on locus ", loc)
      T_stats_loc <- tsum_statistic_1locus(
        strscore_loc = strscore[loc],
        min.quant = min.quant,
        case_control = case_control, trim = trim,
        give.pvalue = give.pvalue, B = B,
        early_stop = early_stop, early_A = early_A, min_stop = early_stop_min
      )
      
      T_stats_list[[loc]] <- T_stats_loc
    }
  }
  T_stats <- rbindlist(T_stats_list, idcol = "locus")
  
  # Prepare output
  outtsum <- exstra_tsum_new_(strscore, tsum = T_stats,
                              correction = correction,
                              alpha = alpha, 
                              args = list(trim = trim, min.quant = min.quant, B = B))
  
  # TODO: remove following lines, maybe
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

# Determine tsum_statistic_1locus assuming input data is for only one locus
# 
# Private function
# 
# @param strscore_loc exstra_strscore object for a single locus (do not give multiple loci).
# @param case_control Set to TRUE if this is a case-control analysis, and FALSE if looking for outliers in all samples.
# @param min.quant Only quantiles above this value are used in constructing the statistic.  
# @param trim Trim this proportion of data points at each quantile 
#             level (rounded up). Must be at least 0 and less than 0.5, but values close
#             to 0.5 may remove all samples and hence result in an error.  
# @param give.pvalue If TRUE, generate the p-value in the simulation.
# @param B Number of simulations in calculating null distributions. The denominator will
#          be B + 1, hence values of B = 10^i - 1 will result in p-values that are 
# @param parallel Use the parallel package when simulating the distribution, creating the
#                 required cluster. 
#                 If cluster is specified then this option makes no difference. 
# @param cluster  A cluster object from the parallel package. Use if you wish to set up 
#                 the cluster yourself or reuse an existing cluster. 
# @param early_stop Simulation may use less replicates when all p-values are large, controlled with early_A.
# @param early_A  Simulations may stop when p.value.sd < early_A * min(p.value). 
#                 Checked approximately when the number of simulations has doubled.
# @param min_stop Minimum number of simulations before stopping.
# @param verbose  If TRUE, report if stopped early, and to what number of simulations.
# 
# @import data.table
# @import stringr
# @import testit
#' @importFrom stats mad median qnorm rnorm
tsum_statistic_1locus <- function(
  strscore_loc,
  case_control = FALSE,
  min.quant = 0,
  trim = 0,
  give.pvalue = FALSE,
  B = 999,
  parallel = FALSE,
  cluster = NULL,
  early_stop = TRUE,
  early_A = 0.25,
  min_stop = 50,
  verbose = TRUE)
{
  qm <- make_quantiles_matrix(strscore_loc, sample = NULL, 
    method = "quantile7")
  
  # trim, min.quant, case_control not default
  # subject is in background by default
  # Use mean and variance of trimmed data (no correction)
  qmt <- qm$y.mat
  # remove lower quantiles
  if(min.quant != 0) {
    qs <- ceiling(ncol(qmt) * (1 - min.quant)) # number columns to keep
    if(qs < 1) {
      qs <- 1
      warning("Number quantile levels used set to 1, trim may be too high for data or you may be too few data points")
    }
    qmt <- qmt[, seq(ncol(qmt) - qs + 1, ncol(qmt))]
  }
  
  if(case_control) {
    qmtest <- qmt[strscore_loc$samples[group == "case" & ! sample %in% qm$low.count, sample], ]
    qmt_bac <- qmt[strscore_loc$samples[group != "case"& ! sample %in% qm$low.count, sample], ]
  } else {
    qmtest <- qmt
    qmt_bac <- qmt
  }
  
  # trim at each quantile
  qmt_bac_untrim <- qmt_bac
  if(trim != 0) {
    qmt_bac <- trim_matrix_(qmt_bac, trim)
  }
  
  # verify qmt is not too small:
  if(nrow(qmt_bac) < 2) {
    stop("Too few samples (", nrow(qmt_bac), ") left after trimming")
  }
  
  # calculate mean, var
  bac_mu <- colMeans(qmt_bac)
  # sum((qm[,1] - bmu[1]) ^ 2 / (nrow(qm) - 1))
  # TODO: may need correction as used in t.test?
  bac_var <- colMeans(qmt_bac ^ 2) - bac_mu ^ 2
  # Made to match R's t.test() functionality:
  bac_stderr <- sqrt(bac_var * (nrow(qmt_bac) + 1) / (nrow(qmt_bac) - 1))
  # bac_s <- sqrt((colMeans(qmt_bac ^ 2) - bac_mu ^ 2) * (nrow(qmt_bac) / (nrow(qmt_bac) - 1)))
  # bac_s <- sqrt(colSums(sweep(qmt_bac, 2, bac_mu) ^ 2) / (nrow(qmt_bac) - 1))
  
  tsums <- qm_tsum_stat_bare_(qmtest, bac_mu, bac_stderr)
  
  if(length(qm$low.count) > 0) {
    tsums_low_count <- rep(NA_real_, length(qm$low.count))
    names(tsums_low_count) <- qm$low.count
    tsums <- c(tsums, tsums_low_count)
  }
  
  #--# P-value generation by simulation #--#
  if(give.pvalue) {
    pvs <- p_value_simulation(tsums = tsums, qmt_bac_untrim = qmt_bac_untrim, 
        case_control = case_control, early_stop = early_stop, min_stop = min_stop,
        B = B, trim = trim, verbose = verbose, early_A = early_A)
    p.value = pvs$p.value
    p.value.sd = pvs$p.value.sd
    B_used = pvs$B_used
  } else {
    p.value = NA_real_ 
    p.value.sd = NA_real_
    B_used = NA_integer_
  }
  
  
  # output data.table directly, so that we can include p-values
  data.table(sample = names(tsums), tsum = tsums, p.value = p.value, p.value.sd = p.value.sd, B = B_used)
}

p_value_simulation <- function(tsums, qmt_bac_untrim, case_control, early_stop, early_A, B, min_stop, trim, verbose){
  tsum_max <- max(tsums, na.rm = TRUE)
  
  # Calculate the mean and se from the median and MAD of the untruncated
  # data instead
  mu_vec <- apply(qmt_bac_untrim, 2, median)
  # NOTE: should the qnorm(3/4) be on the next line? Works better with.
  se_vec <- apply(qmt_bac_untrim, 2, mad) / qnorm(3/4)
  
  N <- nrow(qmt_bac_untrim) # number of samples
  M <- ncol(qmt_bac_untrim) # number of quantiles
  
  # required functions for simulation
  simulate_quantile_matrix <- function() {
    sqm <- t (replicate(N, rnorm(M, mu_vec, se_vec)))
    # as this means the quantiles of a sample are no longer ordered, we sort
    for(i in seq_len(N)) {
      sqm[i, ] <- sort.int(sqm[i,  ], method = "shell")
    }
    sqm
  }
  
  # simulation function for sample in sample in background testing
  sim_tsum_stat_backg <- function() {
    simu <- simulate_quantile_matrix()
    
    # Use same names to make copy-pasta easier
    qmt_bac <- simu 
    
    # trim if required
    if(trim != 0) {
      qmt_bac <- trim_matrix_(qmt_bac, trim)
    }
    
    # calculate mean, var
    bac_mu <- colMeans(qmt_bac)
    bac_var <- colMeans(qmt_bac ^ 2) - bac_mu ^ 2
    # Made to match R's t.test() functionality:
    bac_stderr <- sqrt(bac_var * (nrow(qmt_bac) + 1) / (nrow(qmt_bac) - 1))
    
    tsums <- qm_tsum_stat_bare_(simu, bac_mu, bac_stderr)
    
    tsums
  }
  
  # simulation function for sample in case-control testing
  sim_tsum_stat_cc <- function() {
    qmt_bac <- simulate_quantile_matrix()
    simu_case <- simulate_quantile_matrix()
    
    # trim if required
    if(trim != 0) {
      qmt_bac <- trim_matrix_(qmt_bac, trim)
    }
    
    # calculate mean, var
    bac_mu <- colMeans(qmt_bac)
    bac_var <- colMeans(qmt_bac ^ 2) - bac_mu ^ 2
    # Made to match R's t.test() functionality:
    bac_stderr <- sqrt(bac_var * (nrow(qmt_bac) + 1) / (nrow(qmt_bac) - 1))
    
    tsums <- qm_tsum_stat_bare_(simu_case, bac_mu, bac_stderr)
    
    tsums
  }
  
  # Use the correct simulation function
  if(case_control) {
    sim_tsum_stat <- sim_tsum_stat_cc
  } else {
    sim_tsum_stat <- sim_tsum_stat_backg 
  }
  
  # The simulation
  if(early_stop) {
    if(B <= min_stop) {
      B_part <- B # Otherwise we will have an empty set 
    } else {
      B_part <- floor(B / (2 ^ (seq(ceiling(log2(B/min_stop/2)), 0, -1))) + 2 * .Machine$double.eps)
    }
  } else {
    B_part <- B
  }
  # TODO: explore ways that use far less memory (don't keep raw tsum results)
  #       Is it even worth it? Only for very large B (B*N >> 1e7)
  sim_T <- rep(NA_real_, B*N) # Don't let p-value get to zero
  B_prev <- 0L
  for(B_used in B_part) {
    for(i in seq(B_prev + 1L, B_used, 1L)) {
      sim_T[((i-1)*N + 1) : (i*N) ] <- sim_tsum_stat()
    }
    N_tss <- B_used * N # number tsums in simulation
    if(early_stop & B_used != B) {
      p_smallest <- (sum(sim_T[1L:N_tss] > tsum_max) + 1) / (N_tss + 1)
      p.value.sd <- p_value_sd_(p_smallest, B_used, N)
      # if(p_smallest - p.value.sd * 2 > 0.01) {
      if(verbose && p.value.sd < early_A * p_smallest) {
        message("    Reduced replicates to ", B_used, ".")
        break
      }
    }
    B_prev <- B_used # prepare for next block
  }

  # p-values
  N_out <- length(tsums)
  p.value <- rep(NA_real_, N_out)
  # calculate p-values
  for(i in seq_len(N_out)) {
    p.value[i] <- (sum(sim_T[1L:N_tss] > tsums[i]) + 1) / (N_tss + 1)
  }
  p.value.sd <- p_value_sd_(p.value, B_used, N)
  
  list(p.value = p.value, p.value.sd = p.value.sd, B_used = B_used)
}


qm_tsum_stat_ <- function(qm) {
  
  bmu <- colMeans(qm)
  # sum((qm[,1] - bmu[1]) ^ 2 / (nrow(qm) - 1))
  # TODO: may need correction as used in t.test?
  bvar <- colMeans(qm ^ 2) - bmu ^ 2
  bs <- sqrt(bvar)
  # possibly slower?
  ## bvar <- colSums(sweep(qm, 2, bmu) ^ 2) / nrow(qm)
  # in case a different denominator is required:
  # bvar <- colSums(sweep(qm, 2, bmu) ^ 2) / (nrow(qm) - 1)
  
  list(bmu, bs)
}

#' simple tsum statistic by known 
#' 
#' Assume input has been checked, distribution trimming, and minimum quant chopped.
#' Background mu and variance is given for each quantile.
#' 
#' @param qm Quantile matrix (must be a matrix)
#' @param bmu Background mu vector
#' @param bvar Background variance vector
#' 
#' @keywords internal
qm_tsum_stat_bare_ <- function(
  qm,
  bmu,
  bstde
)
{
  tsums <- rep(NA_real_, nrow(qm))
  names(tsums) <- rownames(qm)
  # check for 0s
  if(any(bstde == 0)) {
    # If there is at least one 0, then spend time correcting this.
    # We are conservative in the correction; this results in a tsum closer to 0.
    stde0 <- TRUE
    stde0v <- bstde == 0
  } else {
    stde0 <- FALSE
  }
  for(i in seq_len(length(tsums))) {
    t_vec <- (qm[i,] - bmu) / bstde
    if(stde0) {
      t_vec[stde0v] <- 0 # conservative correction for 0 variances
    }
    tsums[i] <- mean(t_vec)
  }
  # This is slower when profiled:
  # For one sample: tsum = (mu0 - mu) / S
  # tsums <- rowMeans(sweep(sweep(qm, 2, bmu), 2, bvar, "/"))
  
  tsums
}

#' P-value standard deviation
#' 
#' Accounts for dependencies within each simulation
#'
#' @param p vector of p-values
#' @param B Number of simulations
#' @param N Number of tested samples in each simulation
#' 
#' @keywords internal
p_value_sd_ <- function(p, B, N) {
  Nsim <- B / (1/N + 1/B)
  sqrt(p * ((Nsim + 2)/(Nsim + 1) - p) / Nsim)
}


# Trim a matrix, without preserving sample order
#' @keywords internal
trim_matrix_ <- function(qm, trim = 0) {
  ti <- trim_index_(nrow(qm), trim)
  apply(qm, 2, sort.int, partial = ti)[trim_vector(nrow(qm), trim), ]
}

