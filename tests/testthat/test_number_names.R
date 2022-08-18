# setup
str_score_num <- copy(exstra_wgs_pcr_2["HD"])
str_score_num_0 <- copy(str_score_num)

str_score_num$samples[, sample := str_extract(sample, "(?!0)\\d+")]
str_score_num$data[, sample := str_extract(sample, "(?!0)\\d+")]
setkey(str_score_num$samples, sample)
setkey(str_score_num$data, locus, sample)

str_score_num_0$samples[, sample := str_extract(sample, "\\d+")]
str_score_num_0$data[, sample := str_extract(sample, "\\d+")]
setkey(str_score_num_0$samples, sample)
setkey(str_score_num_0$data, locus, sample)

# the test
test_that("Can handle sample names that look like numbers", {
    expect_s3_class(tsum_test(str_score_num, B = 1), "exstra_tsum")
    expect_s3_class(tsum_test(str_score_num_0, B = 1), "exstra_tsum")
})

