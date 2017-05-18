[![Build Status](https://travis-ci.org/bahlolab/exSTRa.svg?branch=master)](https://travis-ci.org/bahlolab/exSTRa)

# exSTRa
expanded STR algorithm: detecting expansions with paired-end Illumina sequencing data. 



# Structure
At present, the pipeline requires:
- Paired-end Illumina sequencing. Has been tested with some whole-exome sequencing, and whole-genome sequencing with and without PCR in the library preparation step. (Testing on single-end data has not occured, but we presume it would perform poorly due to the affect on alignment.)
- alignment with bowtie2 in local mode (may work for other aligners and settings, not extensively tested due to computational time)
- Sorting, 
- PCR duplicate marking
- Possibly local realignment (need to test)

A database of repeats is required. (Required format to come, but can be flexible and specified, at least in the Perl script.)

Use the Perl scripts and modules to analyse reads in BAM files. This generates STR counts. 

Read above data into R with this package. Provides an OO S3 interface. 
Currently makes extensive use of the data.table package, and understanding its use may help with this package. 
Use the plot() function to draw ECDFs of the data for visual identification. 
Computational methods for mass-screening coming soon... 

# Development

Run `roxygen2::roxygenise()` on code before builds (or commit). Use roxygen2 for documentation. 

