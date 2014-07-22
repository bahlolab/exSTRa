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
  #class(dt) <- c("strdb", class(dt))
  #dt
  structure(list(db = dt, input_type = input_type), class = c("strdb"))
}

strdb_new <- function() {
  # Creates an empty strdb object (may not have any use)
  x <- data.table()
  class(x) <- "strdb"
  x
}

is.strdb <- function(x) inherits(x, "strdb")

strdb_read <- function(file) {
  # Open up a file to load a STR database object
  if (!is.character(file)) stop("file must be character")
  strdatabase <- switch(str_extract(file, "[^.]*$"), 
    xlsx = strdb_xlsx(file),
    strdb_text(file))
  if(is.null(strdatabase)) {
    stop("Could not determine the file type")
  }
  strdatabase
}

strdb_xlsx <- function(file) {
  if (!is.character(file)) stop("file must be character")
  data <- read.xlsx(file, 1)
  assert("xlsx requires Disease column", ! is.null(data$Disease))
  data$disease.symbol <- sub(".*\\((.*)\\).*", "\\1", data$Disease, perl = T)
  strdb(data, "named")
}

strdb_ucsc <- function(file) {
  if (!is.character(file)) stop("file must be character")
  stop("UCSC strdb reading not yet implemented")
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
