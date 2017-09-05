## After running exSTRa_score_analysis.R
context("Test with a real dataset")
# 

expected_text <- "SCA1 (coding CAG) norm: 30 (91bp) , exp: 39 (117bp)"

test_that("text info has changed", {
  expect_equal(loci_text_info(exstra_wgs_pcr_2, "SCA1"), expected_text)
  expect_equal(loci_text_info(exstra_known, "SCA1"), expected_text)
})

