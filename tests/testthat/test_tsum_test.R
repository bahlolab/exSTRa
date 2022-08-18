test_that("tsum_test() may not be giving T statistics", {
  expect_equal(dim(Tsum_stats$stats)[1], 36L)
  expect_equal(length(unique(Tsum_stats$stats$locus)), 2L)
  expect_snapshot(print(tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6", "SCA1", "FRDA")])))
})
