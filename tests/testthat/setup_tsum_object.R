context("Building tsum object")

set.seed(2017)
#exstra_tsum <- tsum_test(exstra_wgs_pcr_2)

Tsum_stats <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], give.pvalue = FALSE)
