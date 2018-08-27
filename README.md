[![Build Status](https://travis-ci.org/bahlolab/exSTRa.svg?branch=master)](https://travis-ci.org/bahlolab/exSTRa)

# exSTRa
expanded STR algorithm: detecting expansions with paired-end Illumina sequencing data. 

This depends on the Perl package 
[Bio::STR::exSTRa](https://github.com/bahlolab/Bio-STR-exSTRa) 
for processing of BAM files. 

# Installation

For easy installation, run from within R:
```
  # install.packages("devtools") # if devtools is not already installed
  devtools::install_github("bahlolab/exSTRa")
```
  

# Using exSTRa
At present, the pipeline requires:
- Paired-end Illumina sequencing. Has been tested with some whole-exome sequencing, and whole-genome sequencing with and without PCR in the library preparation step. (Testing on single-end data has not occured, but we presume it would perform poorly due to the affect on alignment.)
- Alignment with bowtie2 in local mode (may work for other aligners and settings, not extensively tested due to computational time)
- Sorting 
- PCR duplicate marking (recommended)

A database of repeats is required, with known disorder loci included.
An example script to generate a database of all STRs genome wide, or those in genes that are expressed in the brain, is provide in `inst/tools/prepare_exSTRa_input_db.R`.
These input database files can also be [downloaded from FigShare](https://figshare.com/s/0bf679a187d5f3cc2b2c).

Use the Perl scripts and modules from https://github.com/bahlolab/Bio-STR-exSTRa to analyse reads in BAM files. This generates STR counts. 
In the future this functionality may be included within the R exSTRa package. 

This R package provides an OO S3 interface. 
Currently makes extensive use of the data.table package, and understanding its use may help with this package. 

# Examples

Please see the exSTRa vignette for a workable example; 
this can be viewed in R after installation of the package with
```
vignette("exSTRa")
```
or can be downloaded in the repository from [inst/doc/exSTRa.html](inst/doc/exSTRa.html) (in Github, the vignette cannot be viewed directly and should be downloaded).
 
Other datasets in BAM/CRAM format can be analysed in a similar way after processing with the Perl 
[Bio::STR::exSTRa](https://github.com/bahlolab/Bio-STR-exSTRa) package. 

# Citation

Rick M. Tankard, Mark F Bennett, Peter Degorski, Martin B. Delatycki, 
        Paul J. Lockhart, Melanie Bahlo 
         **Detecting tandem repeat expansions in cohorts sequenced with short-read sequencing data**. 
         *bioRxiv* 157792 (2018);
         doi: https://doi.org/10.1101/157792
