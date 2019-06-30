#### Prepare exSTRa database input files
# 
#    Download required input files listed below and change ".../path/to/..." to correct path
#    
#    Script creates annotated input files for running exSTRa on
#    1) all STRs (2-6 bp) in UCSC genome browers simple repeats track
#    2) STRs in genes that are expressed in at least one brain tissue in the GTEx dataset

# ---------------------------------------------------------------------------------------

# Download UCSC genome browser simple repeats track:
# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz
simpleRepeat_file <- ".../path/to/.../simpleRepeat.txt.gz"

# Download GTEx portal median TPM table from https://gtexportal.org/home/datasets (registration required)
GTEx_median_tpm_file <- ".../path/to/.../GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz"

# Specify median TMP value to use as threshold for which genes are considered expressed in brain
brain_median_tpm_thresh <- 1

# Specify miminum and maximum motif size (in base pairs) to search for
min_motif_size <- 2
max_motif_size <- 6

# Download and install ANNOVAR (http://annovar.openbioinformatics.org/)
table_annovar_script <- ".../path/to/.../table_annovar.pl"
humandb_annovar_dir <- ".../path/to/.../humandb"

# ---------------------------------------------------------------------------------------

# Read simpleRepeat table
simpleRepeat <- readr::read_delim(simpleRepeat_file, delim="\t", col_names=FALSE)
colnames(simpleRepeat) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "period", "copyNum", "consensusSize", "perMatch", "perIndel", "score", "A", "C", "G", "T", "entropy", "sequence")
simpleRepeat <- as.data.frame(simpleRepeat, stringsAsFactors=FALSE)

# Filter based on repeat motif size
simpleRepeat <- simpleRepeat[(simpleRepeat$consensusSize >= min_motif_size) & (simpleRepeat$consensusSize <= max_motif_size), ]

# Download GTEx portal median TPM table from https://gtexportal.org/home/datasets:
# GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz
GTEx_median_tpm <- readr::read_delim(GTEx_median_tpm_file, delim="\t", skip=2)

# Get list of genes expressed with median TPM > brain_median_tpm_thresh in at least one brain tissue
brain_tissues <- GTEx_median_tpm[, startsWith(colnames(GTEx_median_tpm), "Brain")]
brain_expressed_genes <- GTEx_median_tpm$Description[rowSums(brain_tissues > brain_median_tpm_thresh) >= 1]


# Prepare data for ANNOVAR
simpleRepeat_avinput <- simpleRepeat[, c("chrom", "chromStart", "chromEnd")]
simpleRepeat_avinput$ref <- "0"
simpleRepeat_avinput$alt <- "0"
write.table(simpleRepeat_avinput, file="simpleRepeat.avinput", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Run ANNOVAR
system(paste0("perl ", table_annovar_script, " simpleRepeat.avinput ", humandb_annovar_dir, " -buildver hg19 -out simpleRepeat_annovar -remove -protocol refGene -operation g -nastring ."))

# Load ANNOVAR input
simpleRepeat_annovar <- read.delim("simpleRepeat_annovar.hg19_multianno.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Convert input to match exSTRa repeat database
simpleRepeat$locus <- paste(simpleRepeat$chrom, simpleRepeat$chromStart, simpleRepeat$sequence, sep="_")
simpleRepeat$long_name <- paste(simpleRepeat$chrom, simpleRepeat$chromStart, simpleRepeat$sequence, sep="_")
simpleRepeat$motif <- simpleRepeat$sequence
simpleRepeat$hg19_start <- simpleRepeat$chromStart + 1
simpleRepeat$hg19_end <- simpleRepeat$chromEnd
simpleRepeat$STR_size_bp <- round(simpleRepeat$copyNum * nchar(simpleRepeat$motif))
simpleRepeat$gene <- simpleRepeat_annovar$Gene.refGene
simpleRepeat$gene_region <- simpleRepeat_annovar$Func.refGene
simpleRepeat$gene_region <- sapply(strsplit(simpleRepeat$gene_region, ";"), function(x) paste(unique(x), collapse=";"))
simpleRepeat$strand <- "+"
simpleRepeat[, c("OMIM", "inheritance", "location", "norm_low", "norm_up", "aff_low", "aff_up", "aff_more", "score_size", "strcat")] <- "NA"

exSTRa_colnames <- c("locus", "long_name", "gene", "location", "gene_region", "motif", "strand", "chrom", "hg19_start", "hg19_end", "copyNum", "perMatch", "perIndel", "STR_size_bp")

simpleRepeat_exSTRa <- simpleRepeat[, exSTRa_colnames]

write.table(simpleRepeat_exSTRa, file="simpleRepeat_repeat_database.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

system("echo '### exSTRa UCSC simpleRepeat Tandem Repeat Finder track database ###' > simpleRepeat_header.txt")
system(paste0("echo '# Last updated ", Sys.Date(), "' >> simpleRepeat_header.txt"))
system("echo '# Requires: exSTRa 0.8' >> simpleRepeat_header.txt")

system("cat simpleRepeat_header.txt simpleRepeat_repeat_database.txt > exSTRa_simpleRepeat.txt")

# Filter genes expressed in brain
simpleRepeat_exSTRa_brain <- simpleRepeat_exSTRa[simpleRepeat_exSTRa$gene %in% brain_expressed_genes, ]

write.table(simpleRepeat_exSTRa_brain, file="simpleRepeat_repeat_database_brain_expressed.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

system("echo '### exSTRa UCSC simpleRepeat Tandem Repeat Finder track database ###' > simpleRepeat_brain_expressed_header.txt")
system("echo '# Filtered to genes expressed in brain (GTEx median tpm > 1 in at least one brain tissue)' >> simpleRepeat_brain_expressed_header.txt")
system(paste0("echo '# Last updated ", Sys.Date(), "' >> simpleRepeat_brain_expressed_header.txt"))
system("echo '# Requires: exSTRa 0.8' >> simpleRepeat_brain_expressed_header.txt")

system("cat simpleRepeat_brain_expressed_header.txt simpleRepeat_repeat_database_brain_expressed.txt > exSTRa_simpleRepeat_brain_expressed.txt")

system("rm simpleRepeat_annovar.hg19_multianno.txt")
system("rm simpleRepeat.avinput")
system("rm simpleRepeat_repeat_database_brain_expressed.txt")
system("rm simpleRepeat_repeat_database.txt")
system("rm simpleRepeat_header.txt")
system("rm simpleRepeat_brain_expressed_header.txt")

