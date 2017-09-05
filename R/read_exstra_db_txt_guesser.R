# Guess the type of data inside a txt file by its header
#' @export
read_exstra_db_txt_guesser <- function(file, ...) {
  # Scan the header
  filefirst <- read.delim(file, nrows = 1, comment.char = "#") 
  
  if("locus" %in% names(filefirst)) {
    # Using a known database format
    read_exstra_db_known(file, ...)
  } else {
    # otherwise default to UCSC style
    read_exstra_db_ucsc(file, ...)
  }
}
