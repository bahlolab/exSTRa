data(exstra_wgs_pcr_2)
image_directory <- tempfile("example_images")
test_that("plot_multi()", {
  expect_null(
    plot_multi(exstra_wgs_pcr_2[c("HD", "SCA6", "FRDA", "SCA1")], dir = image_directory, 
              prefix = "HiSeqXTen_WGS_PCR_2", plot_types = c(1, 2), alpha_case = 0.2))
})
unlink(image_directory)

image_directory_tsum4 <- tempfile("example_images_tsum4")
test_that("plot_multi() should work on exstra_tsum objects (does not at present", {
  expect_error(
    plot_multi(tsum_4, dir = image_directory_tsum4, 
               prefix = "HiSeqXTen_WGS_PCR_2", alpha_case = 0.2))
})
unlink(image_directory_tsum4)
