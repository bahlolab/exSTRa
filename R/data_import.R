# load datasets

#' @include exstra_db.R
#' @include exstra_score.R

#' @export
knownloci_test_db <- exstra_db_read("~/Documents/Research/repeats/disease_repeats/repeat_disorders_2017_04_10.xlsx")

#   str_score <- exstra_score_read (
#     file = "data/HiSeqXTen_WGS_PCR_2.txt", # created by Perl script (TODO: exact name)
#     #database = "data/repeat_disorders.xlsx", # for more control, use object from exstra_db_read() instead
#     database = "../disease_repeats/repeat_disorders_2017_04_26.xlsx",
#     groups.regex = c(case = "^WGSrpt", control = "^WGSrpt_0[24]$"), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
#     filter.low.counts = TRUE
#   )
# 
# str_score