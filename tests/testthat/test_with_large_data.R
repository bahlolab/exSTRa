## After running exSTRa_score_analysis.R
context("Test with a real dataset")
# 
#TODO: make explicit code for test
show(getwd())
wd0 <- getwd()
setwd("../..")
source("examples/exSTRa_score_analysis.R")
setwd(wd0)

expected_text <- "SCA1 (coding CAG) norm: 30 (91bp) , exp: 39 (117bp)"

test_that("text info has changed", {
  expect_equal(loci_text_info(str_score, "SCA1"), expected_text)
  expect_equal(loci_text_info(str_score$db, "SCA1"), expected_text)
})

