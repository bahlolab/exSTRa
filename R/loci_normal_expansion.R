#' @export
loci_normal_exp <- function(x, locus) {
  # Give the reference or normal size of the STR
  if(is.element("strdata", class(x))) {
    x <- x$db
  }
  testit::assert("The class of x must inherit from exstra_db", inherits(x, "exstra_db"))
  locus.in <- locus
  if(x$input_type == "named") {
    x.info <-  x$db[locus.in == locus]
    #TODO: this is wrong
    testit::assert(paste("The locus", locus, "was not found"), dim(x.info)[1] >= 1)
    testit::assert(paste("There were multiple entries for locus", locus), dim(x.info)[1] <= 1)
    rs.len <- with(x.info, nchar(as.character(motif)))
    # normal.copyNum <- with(x.info, floor(copyNum))
    normal.size.bp <- with(x.info, floor(copyNum * rs.len))
    return(c(normal.size.bp, floor(x.info$rn.unst.low * rs.len)))
  } else if (x$input_type == "ucsc") {
    x.info <- x$db[locus.in == locus] 
    #TODO: this is wrong
    testit::assert(paste("The locus", locus, "was not found"), dim(x.info)[1] >= 1)
    testit::assert(paste("There were multiple entries for locus", locus), dim(x.info)[1] <= 1)
    normal.size.bp <- with(x.info, chromEnd - chromStart)
    return(c(normal.size.bp, NA))
  } else {
    stop("Unrecognised input_type in exstra_db. Got ", x$input_type)
  }
}