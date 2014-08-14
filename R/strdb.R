# The strdb class
# Includes information on STRs, whether they are disease-causing or of a general nature

library(data.table)
library(stringr)
library(xlsx)
library(testit)

strdb <- function(dt, input_type = NULL) {
  # Transforms a data.frame or data.table into a strdb object
  if (!is.data.frame(dt)) stop("dt must be data.frame")
  dt <- data.table(dt)
  dt <- dt[!(is.na(chrom) | is.na(chromStart) | is.na(chromEnd))]
  structure(list(db = dt, input_type = input_type), class = c("strdb"))
}

strdb_new <- function() {
  # Creates an empty strdb object (may not have any use)
  x <- data.table()
  class(x) <- "strdb"
  x
}

is.strdb <- function(x) inherits(x, "strdb")

strdb_read <- function(file, ...) {
  # Open up a file to load a STR database object
  if (!is.character(file)) stop("file must be character")
  strdatabase <- switch(str_extract(file, "[^.]*$"), # note str_extract is a stringr function, can be very confusing
    xlsx = strdb_xlsx(file, ...),
    strdb_ucsc(file, ...))
  if(is.null(strdatabase)) {
    stop("Could not determine the file type")
  }
  strdatabase
}

strdb_xlsx <- function(file, ...) {
  if (!is.character(file)) stop("file must be character")
  data <- read.xlsx(file, 1, ...)
  assert("xlsx requires Disease column", ! is.null(data$Disease))
  data <- replace(data, data == "NA", NA)
  data$disease.symbol <- sub(".*\\((.*)\\).*", "\\1", data$Disease, perl = T)
  names(data)[which(names(data) == "hg19.chrom")] <- "chrom"
  names(data)[which(names(data) == "hg19.start.0")] <- "chromStart"
  names(data)[which(names(data) == "hg19.end")] <- "chromEnd"
  strdb(data, "named")
}

strdb_ucsc <- function(file, header = F, ...) {
  if (!is.character(file)) stop("file must be character")
  header.names <- c("X.bin", "chrom", "chromStart", "chromEnd", "name", "period", "copyNum", "consensusSize", "perMatch", "perIndel", "score", "A", "C", "G", "T", "entropy", "sequence")
  col.classes <- c("integer", "factor", "integer", "integer", "factor", 
                   "integer", "numeric", "integer", "integer", "integer", "integer", 
                   "integer", "integer", "integer", "integer", "numeric", "character"
  )
  data <- read.delim(file, header, colClasses = col.classes, ...)
  if(header) {
    assert("When reading UCSC files with a header, the header lines must be: '#bin' (converted to 'X.bin' by R), 'chrom', 'chromStart', 'chromEnd', 'name', 'period', 'copyNum', 'consensusSize', 'perMatch', 'perIndel', 'score', 'A', 'C', 'G', 'T', 'entropy', 'sequence'", 
           identical(names(data), header.names))
  } else {
    assert("When reading UCSC file without a header, must have 17 columns exactly", dim(data)[2] == 17)
    names(data) <- header.names
  }
  assert("The column classes of the UCSC file are not as expected. Are you sure you are supplying a UCSCS file and that the header setting is correct?", prod(sapply(data, class) == col.classes) == 1)
  data$locus <- with(data, paste0(chrom, ":", chromStart + 1, "-", chromEnd, ":", sequence))
  data <- data.table(data)
  data <- data[nchar(sequence) >= 2 & nchar(sequence) <= 6]
  strdb(data, "ucsc")
}

strdb_text <- function(file) {
  if (!is.character(file)) stop("file must be character")
  stop("Text strdb reading not yet implemented")
}

strloci <- function(...) UseMethod("strloci")

strloci.strdb <- function(strdb) {  
  strdatabase$db$disease.symbol
}



Y <- strdb_read("/Users/tankard/Documents/Research/repeats/disease_repeats/repeat_disorders.xlsx")

class(Y)
