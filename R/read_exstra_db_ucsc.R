#' Read simple repeats track that's been downloaded from UCSC Genome Browser.
#' 
#' @param file File path.
#' @param header Explicitly indicate is there is a header in the file. Determined automatically otherwise.
#' @param ... Additional arguments to read.delim()
#' @export
read_exstra_db_ucsc <- function(file, header, ...) {
  if (!is.character(file)) stop("file must be character")
  header.names <- c("X.bin", "chrom", "chromStart", "chromEnd", "name", "period", "copyNum", "consensusSize", "perMatch", "perIndel", "score", "A", "C", "G", "T", "entropy", "sequence")
  col.classes <- c("integer", "factor", "integer", "integer", "factor", 
    "integer", "numeric", "integer", "integer", "integer", "integer", 
    "integer", "integer", "integer", "integer", "numeric", "character"
  )
  
  # Try to work out if this is a table browser downloaded file with a header
  if(missing(header)) {
    filefirst <- read.table(file, nrows = 1, comment.char = "")
    if(filefirst[1, 1] == "#bin") {
      header <- TRUE 
    } else {
      header <- FALSE
    }
  }
  
  data <- read.delim(file, header, colClasses = col.classes, ...)
  # Check if a header line is probably present
  if(header) {
    testit::assert("When reading UCSC files with a header, the header lines must be: '#bin' (converted to 'X.bin' by R), 'chrom', 'chromStart', 'chromEnd', 'name', 'period', 'copyNum', 'consensusSize', 'perMatch', 'perIndel', 'score', 'A', 'C', 'G', 'T', 'entropy', 'sequence'", 
      identical(names(data), header.names))
  } else {
    testit::assert("When reading UCSC file without a header, must have 17 columns exactly", dim(data)[2] == 17)
    names(data) <- header.names
  }
  testit::assert("The column classes of the UCSC file are not as expected. Are you sure you are supplying a UCSCS file and that the header setting is correct?", prod(sapply(data, class) == col.classes) == 1)
  data$locus <- with(data, paste0(chrom, ":", chromStart + 1, "-", chromEnd, ":", sequence))
  data$motif <- data$sequence
  data <- data.table(data)
  data <- data[nchar(sequence) >= 2 & nchar(sequence) <= 6]
  exstra_db_new_(data, "ucsc")
}
