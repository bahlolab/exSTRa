# Bio::STR::exSTRa

Perl module and script for summarising repeat motifs from BAM files.

# Setup 

We plan to make this into a buildable module

To set up, when the working directory is this directory can be achieved with:

    export PERL5LIB="$PWD/exSTRa/perl/lib/:$PERL5LIB"


Dependencies will need to be installed. You can check what cannot yet be found with:

    perl -c main/exSTRa_score.pl

# Script to process BAM/CRAM files

The script `main/exSTRa_score.pl` can assist with this. See `main/run_strexpansion_score.sh` for an example. 

