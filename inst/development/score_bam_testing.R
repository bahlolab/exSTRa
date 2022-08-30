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
x <- score_bam(bam_file)
