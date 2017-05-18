# Bio::STR

Perl module and script for summarising repeat motifs from BAM files.

# Setup 

We plan to make this into a buildable module

To set up, when the working directory is this directory can be achieved with:

    export PERL5LIB="$PWD/exSTRa/perl/lib/:$PERL5LIB"


Dependencies will need to be installed. You can check what cannot yet be found with:

    perl -c STR.pm
    
# Script to process BAM files

The script `main/extract_real_data_read_score.pl` can assist with this. See `main/run_strexpansion_score.sh` for an example. 
