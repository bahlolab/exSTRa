#' @export
loci_normal_exp <- function(x, locus) {
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