test_that("plot_multi()", {
  expect_equal(unlink("example_images", recursive = TRUE), 0) # Remove directory so there isn't a warning
  expect_null(
    plot_multi(exstra_wgs_pcr_2[c("HD", "SCA6", "FRDA", "SCA1")], dir = "example_images", 
              prefix = "HiSeqXTen_WGS_PCR_2", plot_types = c(1, 2), alpha_case = 0.2))
})

test_that("plot_multi() should work on exstra_tsum objects (does not at present", {
  expect_error(
    plot_multi(tsum_4, dir = "example_images_tsum4", 
               prefix = "HiSeqXTen_WGS_PCR_2", alpha_case = 0.2))
})
