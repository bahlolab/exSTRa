context("Testing that the package handles sample names that look like numbers")


str_score <- exstra_wgs_pcr_2
str_score$samples[, sample := str_extract(sample, "(?!0)\\d+")]
str_score$data[, sample := str_extract(sample, "(?!0)\\d+")]
setkey(str_score$samples, sample)
setkey(str_score$data, locus, sample)

test_that("", {
    expect_output(tsum_test(str_score["HD"], B = 10))
})
