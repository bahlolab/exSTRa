# load datasets

#' @include CLASS_exstra_db.R
#' @include CLASS_exstra_score.R
#' @include read_exstra_db.R
#' @include read_score.R
#' @include filter_low_scores.R

#' @export
# TODO: should be importing text file 
exstra_known <- read_exstra_db(system.file("extdata", "repeat_disorders_2017_05_18_external.xlsx", package = "exSTRa"))

#' @export
exstra_wgs_pcr_2 <- read_score (
    system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"),
    database = exstra_known,
    groups.regex = c(case = "^WGSrpt", control = "^WGSrpt_0[24]$"), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
    filter.low.counts = TRUE
  )

