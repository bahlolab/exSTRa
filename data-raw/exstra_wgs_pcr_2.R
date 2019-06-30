# The WGS_PCR_2 data data read in
exstra_wgs_pcr_2 <- read_score (
    system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), # doesn't work before first install
    database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),
    groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
    filter.low.counts = TRUE
  )

wgs_pcr_2_sex <- c(
  "WGSrpt_02" = "female", 
  "WGSrpt_04" = "female", 
  "WGSrpt_05" = "female", 
  "WGSrpt_07" = "female", 
  "WGSrpt_08" = "female", 
  "WGSrpt_09" = "male", 
  "WGSrpt_10" = "female", 
  "WGSrpt_11" = "male", 
  "WGSrpt_12" = "male", 
  "WGSrpt_13" = "female", 
  "WGSrpt_14" = "female", 
  "WGSrpt_15" = "male", 
  "WGSrpt_16" = "male", 
  "WGSrpt_17" = "male", 
  "WGSrpt_18" = "female", 
  "WGSrpt_19" = "male", 
  "WGSrpt_20" = "male", 
  "WGSrpt_21" = "female"
)

exstra_wgs_pcr_2$samples[names(wgs_pcr_2_sex), sex := wgs_pcr_2_sex]

exstra_wgs_pcr_2$samples[, plotname := sub("WGSrpt_", "", sample)]
