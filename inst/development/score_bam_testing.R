# Run score_bam on the Rsamtools example data
library(Rsamtools)
(bam_file <- system.file("extdata", "ex1.bam", package="Rsamtools"))

which <- IRangesList(seq2=IRanges(c(894), c(910)))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)

bam <- scanBam(bam_file, param=param)


names(bam[[1]])

# 
library(GenomicAlignments)

bamp <- readGAlignmentPairs(bam_file, param = param, with.which_label = TRUE)
bampal <- readGAlignmentsList(bam_file, param = param, with.which_label = TRUE)
bamp[[1]]@elementMetadata$seq

# This looks like it has the info we want
bampal[[1]]
# Found the sequences
bampal[[1]]@elementMetadata$seq
bampal[[1]]@elementMetadata$seq[1]
bampal[[1]]@elementMetadata$seq[2]
bamp[[1]]@elementMetadata$seq

# Now how to do regex on that string?
library(Biostrings)
library(stringr)
str_detect(bampal[[1]]@elementMetadata$seq, "CAG")
str_replace_all(bampal[[1]]@elementMetadata$seq, "CAG", "...")
rev(bampal[[1]]@elementMetadata$seq)

X <- "AAGCAGCAGAA"
(x_1 <- str_replace_all(X, "CAG", "..."))
(x_2 <- str_replace_all(X, "AGC", "..."))
(x_3 <- str_replace_all(X, "GCA", "..."))

(y_1 <- str_count(X, "CAG"))
(y_2 <- str_count(X, "AGA"))
(y_3 <- str_count(X, "GCA"))

purrr::map(as.list(bampal), ~ .x@elementMetadata$seq)

count_list <- list()
for(i in seq_along(bampal)) {
  count_list[[i]] <- str_count(bampal[[i]]@elementMetadata$seq, "CAG")
}
unlist(count_list)


# Our test
bamfiles <- c(WGSrpt_12 = "/stornext/Bioinf/data/Bioinformatics/SNPchipdata/MPS_samples/MCRI/AGIP/kinghorn_2017_01_17/repexp_kinghorn_2017_01_17/devel/repexp_2017_01_17_pipeline/bam_recal/WGSrpt_12_bowtie2_recal.bam",
              WGSrpt_18 = "/stornext/Bioinf/data/Bioinformatics/SNPchipdata/MPS_samples/MCRI/AGIP/kinghorn_2017_01_17/repexp_kinghorn_2017_01_17/devel/repexp_2017_01_17_pipeline/bam_recal/WGSrpt_18_bowtie2_recal.bam"
)
x <- score_bam(bamfiles, exstra_known[c("HD", "SCA1")], sample_names = names(bamfiles))

bam <- x$scores[[1]]

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list_loci <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_dt <- purrr::map(list_loci, ~ as.data.table(do.call("DataFrame", .x)))
names(bam_dt) <- exstra_known[c("HD", "SCA1")]$db$locus

rbindlist(bam_dt, idcol = "locus")

dim(bam_df)
