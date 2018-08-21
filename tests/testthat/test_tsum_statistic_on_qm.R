context("Check calculation of tsum statistic")
# 

# Make a small example matrix
qm_small <- make_quantiles_matrix(exstra_wgs_pcr_2["HD", 1:4], 
  method = "quantile7")$y.mat[, c(49, 31, 99)]

bmu <- colMeans(qm_small)
bvar <- colMeans(qm_small ^ 2) - bmu ^ 2

test_that("qm_tsum_stat_bare_ not creating the correct quantile matrix",
  expect_equal(qm_tsum_stat_bare_(qm_small, bmu, bvar),
    c(WGSrpt_02 = 0.138386565958629, WGSrpt_04 = 0.00590495129951718, 
      WGSrpt_05 = 0.179419067436132, WGSrpt_07 = -0.323710584694278), 
    tolerance = 1e-10)
)


# Example with some quantile levels with zero variance
qm_small <- make_quantiles_matrix(exstra_wgs_pcr_2["EPM1A", 2:5], 
  method = "quantile7")$y.mat[, c(49, 31, 60, 99)]

bmu <- colMeans(qm_small)
bvar <- colMeans(qm_small ^ 2) - bmu ^ 2

test_that("qm_tsum_stat_bare_ not creating the correct quantile matrix with some 0 variances",
  expect_equal(qm_tsum_stat_bare_(qm_small, bmu, bvar),
    c(WGSrpt_04 = 0.190251572327043, WGSrpt_05 = -0.0634171907756811, 
      WGSrpt_07 = -0.0634171907756811, WGSrpt_08 = -0.0634171907756811), 
    tolerance = 1e-10)
)