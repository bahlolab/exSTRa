#!/usr/bin/env perl

=head1

within read STR finding from SAM file format. Must be given as STDIN or uncompressed. 


=cut


use 5.014;
use warnings; 
use Data::Dumper;
use Bio::STR::exSTRa::Detection;
use Getopt::Long; 

my $repeat_unit;
GetOptions ("ru=s" => \$repeat_unit)
    or die("Error in command line arguments\n");

my $repeat_unit_rev = scalar reverse $repeat_unit;
$repeat_unit_rev =~ tr/ACGT/TGCA/;

while(<>) {
    chomp;
    if(/^#/) {
        next; # skip SAM header lines
    }
    my @line = split /\t/;
    my $seq = $line[9]; # this is where the sequence should always be
    my $strdetect = Bio::STR::exSTRa::Detection->new (repeat_unit => $repeat_unit, sequence => $seq);
    my $qualloc = $strdetect->qualloc;
    $strdetect = Bio::STR::exSTRa::Detection->new (repeat_unit => $repeat_unit_rev, sequence => $seq);
    my $qualloc_rev = $strdetect->qualloc;
    #printf "%-8s:%s\t%s\n", $qualloc, $seq, $_;
    #printf "%-8s%s\n", $qualloc, $_;
    #say "$qualloc\t$_";
    say "$qualloc\t$qualloc_rev\t$_";
}

