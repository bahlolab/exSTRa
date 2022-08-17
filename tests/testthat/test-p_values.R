test_that("p_values()", {
  expect_equal(nrow(p_values(tsum_4, correction = "samples")), 72)
  expect_equal(nrow(p_values(tsum_4, only.signif = TRUE, correction = "samples")), 8)
  expect_equal(nrow(p_values(tsum_4, only.signif = TRUE, correction = "bf")), 8)
  expect_equal(nrow(p_values(tsum_4, only.signif = TRUE, correction = "loci")), 8)
  expect_equal(nrow(p_values(tsum_4, only.signif = TRUE, correction = "uncorrected")), 10)
  expect_equal(nrow(p_values(tsum_4, only.signif = TRUE, correction = "loci", alpha = 0.4)), 15)
})
