context("Testing rbind_score* functions")

rbinding <- rbind_score_list(list(exstra_wgs_pcr_2[, 5], exstra_wgs_pcr_2[, 10]))

test_that("rbind_score_list() not giving the expected number of outputs", {
  expect_equal(rbinding$samples[, .N], 2L)
  expect_equal(rbinding$data[, .N], 4082L)
  #TODO on $db
})

rbinding <- rbind_score(exstra_wgs_pcr_2[, 5], exstra_wgs_pcr_2[, 10])

test_that("rbind_score() not giving the expected number of outputs", {
  expect_equal(rbinding$samples[, .N], 2L)
  expect_equal(rbinding$data[, .N], 4082L)
  #TODO on $db
})