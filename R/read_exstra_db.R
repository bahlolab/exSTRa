#' Read a file of known STR loci for analysis
#' 
#' Reads a file formatted for exSTRa that contains information of short-tandem repeat
#' loci that may be expanded. 
#' 
#' Files are assumed to be tab-delimited unless the extension is xlsx in which case the 
#' xlsx package is used (if installed).
#' 
#' @param file Path of the file to be read.
#' @param ... Extra arguments to the functions \code{read.delim} (text) or \code{xlsx::read.xlsx} (for xlsx files).
#' 
#' @return An exstra_db object.
#' 
#' @seealso \code{\link{read_score}}
#' 
#' @examples
#' read_exstra_db(system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa")) 
#' 
#' @export
#' @include read_exstra_db_xlsx.R
#' @include read_exstra_db_ucsc.R
#' @include read_exstra_db_known.R
#' @include read_exstra_db_txt_guesser.R
read_exstra_db <- function(file, ...) {
  # Open up a file to load a STR database object
  if (!is.character(file)) stop("file must be character")
  assert("Input file name is an empty string", file != "")
  strdatabase <- switch(str_extract(file, "[^.]*$"), # note str_extract is a stringr function, can be very confusing
    xlsx = read_exstra_db_xlsx(file, ...),
    read_exstra_db_txt_guesser(file, ...))
  if(is.null(strdatabase)) {
    stop("Could not determine the file type")
  }
  strdatabase
}
