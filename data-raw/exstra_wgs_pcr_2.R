# The WGS_PCR_2 data data read in
exstra_wgs_pcr_2 <- read_score (
    system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), # doesn't work before first install
    database = system.file("extdata", "repeat_expansion_disorders", package = "exSTRa"),
    groups.regex = c(case = "^WGSrpt", control = "^WGSrpt_0[24]$"), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
    filter.low.counts = TRUE
  )
