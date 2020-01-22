context("Building tsum object")

set.seed(2017)
four_loci <- c("HD", "SCA6", "SCA1", "FRDA")
tsum_4 <- tsum_test(exstra_wgs_pcr_2[four_loci])

Tsum_stats <- tsum_test(exstra_wgs_pcr_2[c("HD", "SCA6")], give.pvalue = FALSE)
