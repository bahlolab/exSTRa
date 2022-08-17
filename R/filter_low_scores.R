# old name: rsd_filter_lower_than_expected
#' Filter read scores with lower than expected scores
#' 
#' Under the assumption each base in the sequence is uniform and independent.
#' This may help to filter reads not overlapping the STR.
#' 
#' @param strscore exstra_score object.
#' @export
filter_low_scores  <- function(strscore) {
  strscore$db[, unit_length := nchar(as.character(motif))]
  # set score, want to remove scores that are smaller than expected by chance
  strscore$db[, min_score := unit_length / 4 ^ unit_length]
  small_db <- strscore$db[, list(locus, min_score)]
  setkey(small_db, locus)
  # strscore$data <- strscore$data[prop > strscore$db[as.character(locus), min_score]]
  strscore$data <- strscore$data[small_db][prop > min_score][, min_score := NULL]
  setkey(strscore$data, locus, sample)
  #TODO: check bizarre behaviour of data not printing first time here...
  strscore
}


