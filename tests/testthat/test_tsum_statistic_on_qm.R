context("Check calculation of tsum statistic")
# 

# Make a small example matrix
qm_small <- make_quantiles_matrix(exstra_wgs_pcr_2["HD", 1:4], 
  method = "quantile7")$y.mat[, c(49, 31, 99)]

bmu <- colMeans(qm_small)
bvar <- colMeans(qm_small ^ 2) - bmu ^ 2

test_that("qm_tsum_stat_bare_ not creating the correct quantile matrix",
  expect_equal(qm_tsum_stat_bare_(qm, bmu, bvar),
    c(WGSrpt_02 = -0.15344934400607, WGSrpt_04 = -0.464190506848193, 
      WGSrpt_05 = 0.617639850854263), 
    tolerance = 1e-10)
)