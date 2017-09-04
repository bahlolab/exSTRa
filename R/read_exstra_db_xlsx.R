read_exstra_db_xlsx <- function(file, ...) {
  if (!is.character(file)) stop("file must be character")
  data <- xlsx::read.xlsx(file, 1, stringsAsFactors = FALSE, ...)
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
  names(data)[which(names(data) == "Repeat.sequence")] <- "motif"
  
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
