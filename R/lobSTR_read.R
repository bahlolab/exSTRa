# Functions for reading lobSTR VCF file
library(VariantAnnotation)

lobSTR_read <- function(
    file,
    strcounts
  )
{
  
  
  
}


vcfFile = "/Volumes/tankard/projects/research/STRs/repeat_expansion_rediscovery/lobSTR/str_calls/rediscovery.sort.vcf.gz"
tabixfile <- TabixFile(vcfFile)
scanVcfHeader(tabixfile)

singlestr <- strdatabase$db[Location.of.repeat.within.gene == "coding"][1,]
startLoc <- as.numeric(as.character(singlestr$chromStart))
endLoc <- as.numeric(as.character(singlestr$chromEnd))
param = ScanVcfParam(
  info=c("END", "REF", "RPA", "MOTIF"),
  geno="GT",
  which=GRanges(as.character(singlestr$chrom), 
    IRanges(startLoc, endLoc))
)

STR1 <- readVcf(vcfFile, "hg19", param=param)
STR1 <- STR1[info(STR1)$END == endLoc]
# TODO: check that an STR has been found, and if not maybe report on number of potential STRs (and also report number of total STRs at location found)
# TODO: check that a single STR has been found
STR1 <- STR1[1,]
geno(STR1)
rowData(STR1)
info(STR1)
# TODO: must check we only get a single STR and it is the one we want, as well as filter out others
# info(STR1) holds the key to getting what we want

# allele lengths, with an index on the vector of +1 due to no zero index in R
alleles.ixplus1 <- c(info(STR1)$REF[[1]], info(STR1)$RPA[[1]])
  
geno(STR1)$GT
# TODO: translate to repeat lengths
genotypes <- strsplit(geno(STR1)$GT, "/")
for(g in seq(1, length(genotypes), 1)) {
  this.g <- genotypes[[g]]
  if(identical(this.g, ".")) {
    genotypes[[g]] <- c(NA, NA)
  } else {
    genotypes[[g]] <- alleles.ixplus1[as.numeric(this.g) + 1]
  }
}


