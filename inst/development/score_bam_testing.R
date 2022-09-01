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
cluster25 <- snow::makeCluster(25)
library(fs)
bamfiles <- dir_ls("/stornext/Bioinf/data/Bioinformatics/SNPchipdata/MPS_samples/MCRI/AGIP/kinghorn_2017_01_17/repexp_kinghorn_2017_01_17/devel/repexp_2017_01_17_pipeline/bam_recal/",
       glob = "*.bam")
names(bamfiles) <- str_extract(bamfiles, "(?<=bam_recal/).+(?=_bowtie2_recal)")

x_overlap <- score_bam(bamfiles, exstra_known, sample_names = names(bamfiles), groups.regex = c(case = "^WGSrpt", control = ""), 
               verbosity = 2, filter.low.counts = TRUE, cluster = cluster25)
x_count <- score_bam(bamfiles, exstra_known, sample_names = names(bamfiles), groups.regex = c(case = "^WGSrpt", control = ""), 
                       verbosity = 2, filter.low.counts = TRUE, method = "count", cluster = cluster25)

tsxo <- tsum_test(x_overlap, correction = "samples", cluster = cluster25)
tsxc <- tsum_test(x_count, B = 999, correction = "samples", cluster = cluster25)
tswgs2 <- tsum_test(exstra_wgs_pcr_2)

par(mfrow = c(3, 1))
plot(tsxc["HD"])
plot(tsxo["HD"])
plot(tswgs2["HD"])

p_values(tsxo, only.signif = TRUE)
p_values(tsxc, only.signif = TRUE)
p_values(tswgs2, only.signif = TRUE)

par(mfrow = c(1, 1))
plot(tsxc)

sbf <- scanBamFlag(
  isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
  isNotPassingQualityControls = FALSE, isDuplicate = FALSE
)

Y <- score_bam_1(bamfiles[1], exstra_known, sample_names = names(bamfiles)[1], scan_bam_flag = sbf)


bamheader <- scanBamHeader(bamfiles[[1]])


x_count1 <- score_bam(bamfiles, exstra_known["SCA1"], sample_name_origin = "basename", 
                      sample_name_remove = "_bowtie2_recal",
                      groups.regex = c(case = "^WGSrpt", control = ""), 
                     verbosity = 2, filter.low.counts = TRUE, method = "count")


# Timing parallel implementation
start.time <- Sys.time()
x_count_single <- score_bam(bamfiles, exstra_known, sample_name_origin = "basename", 
                      sample_name_remove = "_bowtie2_recal",
                      groups.regex = c(case = "^WGSrpt", control = ""), 
                      verbosity = 2, filter.low.counts = TRUE, method = "count")
end.time <- Sys.time()
(time.taken_single <- end.time - start.time)

start.time <- Sys.time()
x_count_parallel <- score_bam(bamfiles, exstra_known, sample_name_origin = "basename", 
                      sample_name_remove = "_bowtie2_recal",
                      groups.regex = c(case = "^WGSrpt", control = ""), 
                      verbosity = 2, filter.low.counts = TRUE, method = "count", parallel = TRUE)
end.time <- Sys.time()
(time.taken_parallel <- end.time - start.time)


start.time <- Sys.time()
x_count_parallel <- score_bam(bamfiles, exstra_known, sample_name_origin = "basename", 
                              sample_name_remove = "_bowtie2_recal.bam$",
                              groups.regex = c(case = "^WGSrpt", control = ""), 
                              verbosity = 2, filter.low.counts = TRUE, method = "count", 
                              parallel = TRUE, cluster = cluster25)
end.time <- Sys.time()
(time.taken_parallel_premade_cluster <- end.time - start.time)


start.time <- Sys.time()
tst_count_single <- tsum_test(x_count_parallel, B = 9999, correction = "samples")
end.time <- Sys.time()
(time.taken_tst_count_single <- end.time - start.time)

start.time <- Sys.time()
tst_count_parallel <- tsum_test(x_count_parallel, B = 9999, cluster = cluster25, correction = "samples")
end.time <- Sys.time()
(time.taken_tst_count_parallel_premade_cluster <- end.time - start.time)

start.time <- Sys.time()
tst_count_snow <- tsum_test_snow(x_count_parallel, B = 9999, cluster = cluster25, correction = "samples")
end.time <- Sys.time()
(time.taken_tst_count_snow_premade_cluster <- end.time - start.time)


# getting current score data
# Read score data and file with loci information
str_score_wgs_pcr_2 <- read_score (
  file = "/stornext/Bioinf/data/lab_bahlo/users/tankard.r/exstra/wgs_pcr_2/wgs_pcr_2_exSTRa_scores.txt", 
  database = system.file("extdata", "repeat_expansion_disorders_hg19.txt", package = "exSTRa"),  # for greater control, use object from read_exstra_db() instead
  groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), # the group is the first regular expression (regex) to match
  filter.low.counts = TRUE
)

tsxo <- tsum_test(x_overlap, correction = "samples", cluster = cluster25)
tsxc <- tsum_test(x_count, B = 999, correction = "samples", cluster = cluster25)
ts_wgs_pcr_2 <- tsum_test(str_score_wgs_pcr_2, correction = "samples", cluster = cluster25)

par(mfrow = c(3, 1))
plot(tsxc["HD"], main = "HD, scoring by count")
plot(tsxo["HD"], main = "HD, scoring by overlap (R)")
plot(ts_wgs_pcr_2["HD"], main = "HD, scoring by overlap (Perl)")

p_values(tsxc, only.signif = TRUE)

tsxc_bf <- copy(tsxc)
p_values(tsxc_bf, modify = TRUE)
tsxo_bf <- copy(tsxo)
p_values(tsxo_bf, modify = TRUE)
ts_wgs_pcr_2_bf <- copy(ts_wgs_pcr_2)
p_values(ts_wgs_pcr_2_bf, modify = TRUE)
par(mfrow = c(3, 1))
plot(tsxc["HD"], main = "HD, scoring by count")
plot(tsxo["HD"], main = "HD, scoring by overlap (R)")
plot(ts_wgs_pcr_2["HD"], main = "HD, scoring by overlap (Perl)")
