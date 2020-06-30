#' @import data.table
#' @import testit
#' @import checkmate
read_exstra_db_known <- function(file, ...) {
  if (!is.character(file)) stop("file must be character")
  data <- read.delim(file, stringsAsFactors = FALSE, comment.char = "#", ...)
  assert("Known requires Disease or locus column", ! is.null(data$Disease) || ! is.null(data$locus))
  if(is.null(data$Disease)) {
    assert("Known requires Disease or locus column", ! is.null(data$locus))
    data$Disease <- data$locus
  }
  data <- replace(data, data == "NA", NA)
  data$locus <- sub(".*\\((.*)\\).*", "\\1", data$Disease, perl = T)
  
  # Match the first suitable column
  names(data)[assert_int(grep("chr(om)?$", names(data), TRUE))] <- "chrom"
  names(data)[assert_int(grep("start(\\.0)?$", names(data), TRUE))] <- "chromStart"
  names(data)[assert_int(grep("end$", names(data), TRUE))] <- "chromEnd"
  
  # give more verbose repeat number information
  
  exstra_db_new_(data, "named")
}
