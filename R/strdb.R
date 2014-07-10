# The strdb class
# Includes information on STRs, whether they are disease-causing or of a general nature

library(data.table)
library(stringr)
library(xlsx)

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
  strdb(read.xlsx(file, 1))
}

strdb_ucsc <- function(file) {
  if (!is.character(file)) stop("file must be character")
  stop("UCSC strdb reading not yet implemented")
}

strdb_text <- function(file) {
  if (!is.character(file)) stop("file must be character")
  stop("Text strdb reading not yet implemented")
}

strdb_new <- function() {
  x <- data.table()
  class(x) <- "strdb"
  x
}

strdb <- function(dt) {
  if (!is.data.frame(dt)) stop("dt must be data.frame")
  dt <- data.table(dt)
  class(dt) <- c("strdb", class(dt))
  dt
  #structure(list(db = dt), class = c("strdb"))
}

is.strdb <- function(x) inherits(x, "strdb")


Y <- strdb_read("/Users/tankard/Documents/Research/repeats/disease_repeats/repeat_disorders.xlsx")

class(Y)
