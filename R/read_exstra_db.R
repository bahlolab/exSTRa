#' @export
#' @include read_exstra_db_xlsx.R
#' @include read_exstra_db_ucsc.R
read_exstra_db <- function(file, ...) {
  # Open up a file to load a STR database object
  if (!is.character(file)) stop("file must be character")
  strdatabase <- switch(str_extract(file, "[^.]*$"), # note str_extract is a stringr function, can be very confusing
    xlsx = read_exstra_db_xlsx(file, ...),
    read_exstra_db_ucsc(file, ...))
  if(is.null(strdatabase)) {
    stop("Could not determine the file type")
  }
  strdatabase
}
