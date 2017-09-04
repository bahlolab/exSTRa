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
  names(data)[which(names(data) == "hg19.chrom" | names(data) == "hg19_chr")] <- "chrom"
  names(data)[which(names(data) == "hg19.start.0" | names(data) == "repeat.start" | names(data) == "hg19_start")] <- "chromStart"
  names(data)[which(names(data) == "hg19.end" | names(data) == "repeat.end" | names(data) == "hg19_end")] <- "chromEnd"
  
  # give more verbose repeat number information
  
  exstra_db_new_(data, "named")
}
