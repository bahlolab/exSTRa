# Preparing BAM for exSTRa

*19 May 2017 12:07 AM AEST*

This is a guide on how we have analysed our FastQ data in order to be ready for the Perl script of exSTRa. 
Sequencing data should be Illumina paired-end. 

At present we have only tested our software with Bowtie2 (2.2.9) with `--very-sensitive-local --maxins 1000`.
It is possible our software will work satisfactorily with other aligners and settings; we plan to test some common possibilities in the near future including bwa-mem, novoalign, and bowtie2 with default settings. 

See the end of this page for version numbers.

## Alignment 

### Setting readgroups

Set the variables with $ as is expected for BAM readgroups. 
Our Perl script will use the SM (sample) tag in the output as the sample name. 
Our Perl script only uses --rg-id (ID) and SM. 
Other tags may affect other tools, so set according to BAM specification. 
If you have more than 

    bowtie2_readgroups="--rg-id $RGID --rg PL:ILLUMINA --rg PU:$PlatformUnit --rg LB:$Library --rg DS:$Description --rg SM:$sample --rg CN:$centre"


### bowtie2 alignment

Create a bowtie2 index for a hg19 reference sequence in the order expected for GATK. 
(see [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle), file `bundle/2.3/hg19/ucsc.hg19.fasta` (later versions should be identical (TODO verify))).
Please verify `--phred33` is appropriate for your data. 
We run bowtie2 on each pair of FastQ files.

    # 
    reference=hg19.fa # reference FASTA path
    bowtie2_ref=hg19.fa # bowtie2 index prefix path
    fastq1= # FastQ read 1 path
    fastq2= # FastQ read 2 path
    align_cores=40 # alignment cores
    outbam=aligned_1.bam # output BAM

    bowtie2 \
        -x $bowtie2_ref \
        -1 ${fastq1} \
        -2 ${fastq2} \
        --phred33 \
        --threads $align_cores \
        --very-sensitive-local \
        --maxins 1000 \
        $bowtie2_readgroups \
        2> ${bamalign}.bowtie2.stderr.txt \
        | samtools view -SbT ${reference} - \
        > ${bamalign}  

## Sorting and merge

The choice of sorting software should be mostly inconsequential. 
Sorting was performed for each sample, merging separate alignments of the same sample.
Our pipeline has used Novosort, as follows:

    sorting_cores= # number of cores for sorting. Without a licence this is limited to 1
    tmpdir= # temporary directory, should be much larger than the total size of input BAMs
    
    # make sure the RAM is appropriate for your system, 
    # makes the biggest difference if this can hold the uncompressed BAM all in memory
    novosort \
    --threads "$sorting_cores" \
    --tmpdir "${TempDir}" \
    --index \
    --ram 100G \
    --output output.sort.bam \
    aligned_1.bam [aligned_2.bam] [...]

### Duplicate marking

Duplicate removal

    PICARDPATH=picard.jar # path to picard.jar
    # Assuming 16GB+ memory on your machine
    java -jar -Xmx10g -XX:ParallelGCThreads=1 $PICARDPATH MarkDuplicates


# Versions of software tested

All software was run on Linux servers. 

[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): 2.2.9

[samtools](http://www.htslib.org/doc/samtools.html) (untracked)

[novosort](http://www.novocraft.com/products/novosort/) (3.04.04)

[Picard](https://broadinstitute.github.io/picard/) 2.1.1(6a5237c0f295ddce209ee3a3a5b83a3779408b1b_1457101272) 

java

    java version "1.8.0_40"
    Java(TM) SE Runtime Environment (build 1.8.0_40-b26)
    Java HotSpot(TM) 64-Bit Server VM (build 25.40-b25, mixed mode)
