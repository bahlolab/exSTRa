# load datasets

#' @include CLASS_exstra_db.R
#' @include CLASS_exstra_score.R

#' @export
# TODO: should be importing text file 
exstra_known <- exstra_db_read("~/Documents/Research/repeats/disease_repeats/repeat_disorders_2017_04_10.xlsx")

#' @export
exstra_wgs_pcr_2 <- exstra_score_read (
    system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"),
    database = exstra_known,
    groups.regex = c(case = "^WGSrpt", control = "^WGSrpt_0[24]$"), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
    filter.low.counts = TRUE
  )

