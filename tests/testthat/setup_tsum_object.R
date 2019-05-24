context("Building tsum object")

set.seed(2017)
tsum_4 <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6", "SCA1", "FRDA")])

Tsum_stats <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], give.pvalue = FALSE)
