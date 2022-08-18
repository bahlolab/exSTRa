test_that("tsum_plots", {
  vdiffr::expect_doppelganger("SCA3 ECDF", fig = function() plot(exstra_wgs_pcr_2["SCA3"]))
})