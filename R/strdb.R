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
  assert("xlsx requires Disease or locus column", ! is.null(data$Disease) || ! is.null(data$locus))
  if(is.null(data$Disease)) {
    assert("xlsx requires Disease or locus column", ! is.null(data$locus))
    data$Disease <- data$locus
  }
  data <- replace(data, data == "NA", NA)
  data$disease.symbol <- sub(".*\\((.*)\\).*", "\\1", data$Disease, perl = T)
  names(data)[which(names(data) == "hg19.chrom" | names(data) == "hg19_chr")] <- "chrom"
  names(data)[which(names(data) == "hg19.start.0" | names(data) == "repeat.start")] <- "chromStart"
  names(data)[which(names(data) == "hg19.end" | names(data) == "repeat.end")] <- "chromEnd"
  
  # give more verbose repeat number information
  data$rn.stab.low <- as.numeric(NA)
  data$rn.stab.hig <- as.numeric(NA)
  data$rn.unst.low <- as.numeric(NA) 
  data$rn.unst.hig <- as.numeric(NA) 
  data$rn.unst.nonmax <- F
  
  for(i in 1:dim(data)[1]) {
    data[i, c("rn.stab.low", "rn.stab.hig")] <- as.numeric(strsplit(as.character(data[i, "Stable.repeat.number"]), "-")[[1]])
    if(grepl("\\+", as.character(data[i, "Unstable.repeat.number"]))) {
      data[i, "rn.unst.nonmax"] <- T
    }
    if(grepl("-", as.character(data[i, "Unstable.repeat.number"]))) {
      data[i, c("rn.unst.low", "rn.unst.hig")] <- as.numeric(strsplit(sub("\\+", "", as.character(data[i, "Unstable.repeat.number"])), "-")[[1]])
    } else {
      data[i, c("rn.unst.low", "rn.unst.hig")] <- rep(as.numeric(sub("\\+", "", as.character(data[i, "Unstable.repeat.number"]))), 2)
    }
  }
  
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
  loci <- strdb$db$disease.symbol
  if(is.null(loci)) {
    loci <- strdb$db$locus
  }
  assert("Could not identify the str loci", !is.null(loci))
  loci
}

strloci_text_info <- function(strdata, locus) {
  # gives text info for the locus, usually used in plot titles
  
}


Y <- strdb_read("/Users/tankard/Documents/Research/repeats/disease_repeats/repeat_disorders.xlsx")

class(Y)
