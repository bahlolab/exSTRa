library(data.table)
library(exSTRa)
library(glue)

# Write bed file of reads to extract
known_extraction_bed <- exstra_known$db[, .(chrom, start = chromStart - 1001, end = chromEnd + 1000)]
readr::write_tsv(known_extraction_bed, path = "inst/extdata/known_bam_extract.bed")

# Extract reads using external tool samtools
bamfile <- "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam"
system(glue("bash samtools view -L inst/extdata/known_bam_extract.bed {bamfile}"))

# Also do this with:
# ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam
# ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam

# Filter reads that do not have both pairs (preventing warnings)


