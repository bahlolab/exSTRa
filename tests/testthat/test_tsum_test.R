data(exstra_wgs_pcr_2)
test_that("tsum_test()", {
  expect_equal(dim(Tsum_stats$stats)[1], 36L)
  expect_equal(length(unique(Tsum_stats$stats$locus)), 2L)
  expect_snapshot(print(tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6", "SCA1", "FRDA")])))
})

test_that("tsum_plots", {
  vdiffr::expect_doppelganger("HD ECDF", fig = function() plot(tsum_4["HD"]))
})
