# strexpansion
str expansion R package


# Package name

There is potential to change the package name. There is a brainstorm of ideas in 
[alternative_names.txt](alternative_names.txt)


# Structure
At present, the pipeline requires:
- alignment with bowtie2 in local mode (may work for other aligners and settings, not extensively tested due to computational time)
- Sorting, 
- PCR duplicate marking
- Possibly local realignment (need to test)

A database of repeats is required. (Required format to come, but can be flexible and specified, at least in Perl)

Then use Perl scripts and modules to analyse reads in BAM files. This generates STR counts. 

Read above data into R with this package. Provides an OO S3 interface. 
Currently makes extensive use of the data.table package, and understanding its use may help with this package. 
Use the plot() function to draw ECDFs of the data for visual identification. 
Computational methods for mass-screening coming soon... 
