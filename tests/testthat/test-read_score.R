test_that("read_score()", {
  expect_true(
    inherits(
      read_score(
        file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
        database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),  # for greater control, use object from read_exstra_db() instead
        groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), # the group is the first regular expression (regex) to match
        filter.low.counts = TRUE
      ),
    "exstra_score")
  )
})
