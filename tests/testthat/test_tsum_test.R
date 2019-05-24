context("Testing tsum_test() functions")

test_that("tsum_test() may not be giving T statistics", {
  expect_equal(dim(Tsum_stats$stats)[1], 36L)
  expect_equal(length(unique(Tsum_stats$stats$locus)), 2L)
})
