# The exstra_db class
# Includes information on STRs, whether they are disease-causing or of a general nature

# check if the object is of this class
#' @import data.table
#' @import stringr
#' @import xlsx
#' @import testit
#' 
#' @export
is.exstra_db <- function(x) inherits(x, "exstra_db")

# Create a new object of this class (not for the user)
#
exstra_db_new_ <- function(strd, input_type = NULL) {
  # Transforms a data.frame or data.table into a exstra_db object
  if (!is.data.frame(strd)) stop("strd must be data.frame")
  strd <- data.table(strd)
  strd <- strd[!(is.na(chrom) | is.na(chromStart) | is.na(chromEnd))]
  strd$input_order <- seq(1, dim(strd)[1], 1)
  # TODO this should always just be locus, hack away!
  if(!is.null(strd$disease.symbol)) {
    #setkey(strd, "disease.symbol")
    setnames(strd, "disease.symbol", "locus")
  } 
  setkey(strd, "locus")
  structure(list(db = strd, input_type = input_type), class = c("exstra_db"))
}

#' @export
exstra_db_read <- function(file, ...) {
  # Open up a file to load a STR database object
  if (!is.character(file)) stop("file must be character")
  strdatabase <- switch(str_extract(file, "[^.]*$"), # note str_extract is a stringr function, can be very confusing
    xlsx = exstra_db_xlsx(file, ...),
    exstra_db_ucsc(file, ...))
  if(is.null(strdatabase)) {
    stop("Could not determine the file type")
  }
  strdatabase
}

exstra_db_xlsx <- function(file, ...) {
  if (!is.character(file)) stop("file must be character")
  data <- read.xlsx(file, 1, stringsAsFactors = FALSE, ...)
  assert("xlsx requires Disease or locus column", ! is.null(data$Disease) || ! is.null(data$locus))
  if(is.null(data$Disease)) {
    assert("xlsx requires Disease or locus column", ! is.null(data$locus))
    data$Disease <- data$locus
  }
  data <- replace(data, data == "NA", NA)
  data$locus <- sub(".*\\((.*)\\).*", "\\1", data$Disease, perl = T)
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
  
  exstra_db_new_(data, "named")
}

#' @export
exstra_db_ucsc <- function(file, header = F, ...) {
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
  data$Repeat.sequence <- data$sequence
  data <- data.table(data)
  data <- data[nchar(sequence) >= 2 & nchar(sequence) <= 6]
  exstra_db_new_(data, "ucsc")
}

#' @export
print.exstra_db <- function(x, ...) {
  cat(class(x)[1], " object with ", dim(x$db)[1], " loci ($db) of type ",  x$input_type, "\n",
    sep = "")
}

#' @export
exstra_db_text <- function(file) {
  if (!is.character(file)) stop("file must be character")
  stop("Text exstra_db reading not yet implemented")
}

#' @export
loci.exstra_db <- function(exstra_db) {
  # Give the loci names
  loci <- exstra_db$db[order(input_order), locus]
  assert("Could not identify the str loci", !is.null(loci))
  loci
}

#' @export
loci_text_info.exstra_db <- function(x, locus) {
  # gives text info for the locus, usually used in plot titles
  # TODO: modify this:
  assert("The class of x must be exstra_db", is.exstra_db(x))
  locus.in <- locus
  if(x$input_type == "named") {
    x.info <- x$db[locus.in == locus]
    #TODO: this is wrong
    assert(paste("The locus", locus, "was not found"), dim(x.info)[1] >= 1)
    assert(paste("There were multiple entries for locus", locus), dim(x.info)[1] <= 1)
    rs.len <- with(x.info, nchar(as.character(Repeat.sequence)))
    normal.copyNum <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum), floor(read_detect_size / rs.len)))
    normal.size.bp <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum * rs.len), read_detect_size))
    return(with(x.info,  
      paste0(locus, " (", 
        Location.of.repeat.within.gene, " ", Repeat.sequence, ") norm: ", normal.copyNum, 
        " (", normal.size.bp, "bp) , exp: ", rn.unst.low, " (", 
        floor(rn.unst.low * rs.len), "bp)"))
    )
  } else if (x$input_type == "ucsc") {
    x.info <- x$db[locus.in == locus] 
    #TODO: this is wrong
    assert(paste("The locus", locus, "was not found"), dim(x.info)[1] >= 1)
    assert(paste("There were multiple entries for locus", locus), dim(x.info)[1] <= 1)
    rs.len <- with(x.info, nchar(as.character(Repeat.sequence)))
    #normal.copyNum <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum), floor(read_detect_size / rs.len)))
    #normal.size.bp <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum * rs.len), read_detect_size))
    return(with(x.info,  
      paste0(locus))
    )
  } else {
    stop("Unrecognised input_type in exstra_db. Got ", x$input_type)
  }
}

#' @export
strloci_normal_exp <- function(x, locus) {
  # Give the reference or normal size of the STR
  if(is.element("strdata", class(x))) {
    x <- x$db
  }
  assert("The class of x must be exstra_db or strdata", class(x) == "exstra_db")
  locus.in <- locus
  if(x$input_type == "named") {
    x.info <-  x$db[locus.in == locus]
    #TODO: this is wrong
    assert(paste("The locus", locus, "was not found"), dim(x.info)[1] >= 1)
    assert(paste("There were multiple entries for locus", locus), dim(x.info)[1] <= 1)
    rs.len <- with(x.info, nchar(as.character(Repeat.sequence)))
    normal.copyNum <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum), floor(read_detect_size / rs.len)))
    normal.size.bp <- with(x.info, ifelse(is.null(read_detect_size), floor(copyNum * rs.len), read_detect_size))
    return(c(normal.size.bp, floor(x.info$rn.unst.low * rs.len)))
  } else if (x$input_type == "ucsc") {
    x.info <- x$db[locus.in == locus] 
    #TODO: this is wrong
    assert(paste("The locus", locus, "was not found"), dim(x.info)[1] >= 1)
    assert(paste("There were multiple entries for locus", locus), dim(x.info)[1] <= 1)
    stop("Not programmed for UCSC input")
    return(NULL
    )
  } else {
    stop("Unrecognised input_type in exstra_db. Got ", x$input_type)
  }
}

#' @export
strloci_normal <- function(x, locus) {
  strloci_normal_exp (x, locus)[1]
}

#' @export
strloci_minexp <- function(x, locus) {
  # Give the minimum expanded STR in bp
  strloci_normal_exp (x, locus)[2]
}

#' @export
`[.exstra_db` <- function(x, fil) {
  assert("locus not the key of x$db", key(x$db)[1] == "locus")
  x$db <- x$db[eval(substitute(fil))]
  x
}

# I think the following was code that was left over from another time
#Y <- exstra_db_read("/Users/tankard/Documents/Research/repeats/disease_repeats/repeat_disorders.xlsx")
#class(Y)


# TODO method for seqnames(exstra_db)