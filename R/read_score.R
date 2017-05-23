#' @export
read_score <- function(file, database, groups.regex = NULL, groups.samples = NULL, filter.low.counts = TRUE) {
  # Load the STR counts
  # Groups should be named null, control and case
  if(is.character(database)) {
    # as database is presumbly a file, try to read from it
    database <- read_exstra_db(database)
  }
  out <- strs_read_(file, database, groups.regex, groups.samples, this.class = "exstra_score")
  # TODO: checks for this input
  out$data$prop <- with(out$data, rep / mlength)
  strscore <- exstra_score_new_(out$data, out$db)
  if(filter.low.counts) {
    # Filter low counts, assumed wanted by default
    strscore <- filter_low_scores(strscore)
  }
  return(strscore)
}
