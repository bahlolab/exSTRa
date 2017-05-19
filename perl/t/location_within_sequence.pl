#!/usr/bin/env perl

=head1

within read STR finding

Some ideas to find STRs within reads. 

=cut


use 5.014;
use warnings; 
use Data::Dumper;
use Bio::STR::exSTRa::Detection;

my %headcols;
my ($tests_passed, $n_tests) = (0, 0);
my $test_qualloc = '';
my $qualloc;
while (<>) {
    chomp;
    if(/^#END#/) {
        last; # don't bother with any more columns
    }
    s/"//g; # assuming we aren't really using these quotes
    my @line = split /\t/;
    if($. == 1) {
        @headcols{@line} = (0..$#line);
        if (defined ($headcols{qual_loc})) {
            $test_qualloc = 1;
        }
        next;
    }
    my ($rep_unit, $seq, $starts, $ends, $lengths) = @line[@headcols{qw(repeat_unit sequence starts ends lengths)}];
    if(!defined($rep_unit) || !defined($seq) || !defined($starts) || !defined($ends) || !defined($lengths)) {
        die "There does not appear to be all of the columns 'repeat_unit', 'sequence', 'starts', 'ends' and 'lengths' in the input\n";
    }
    my $strdetect = Bio::STR::exSTRa::Detection->new (repeat_unit => $rep_unit, sequence => $seq);
    $strdetect->verify_repeat_found;
    say $line[$headcols{Disease}];
    say $strdetect->{sequence};
    my $location_info = $strdetect->location;
    for my $key (keys %$location_info) {
        $location_info->{$key} = join(",", @{$location_info->{$key}});
    }
    if($location_info->{starts} eq $starts && $location_info->{ends} eq $ends && $location_info->{lengths} eq $lengths) {
        $tests_passed++;
        say "Location matches";
    } else {
        say "Location does not match:";
        say "starts: file: $starts; sequence: $location_info->{starts}";
        say "ends: file: $ends; sequence: $location_info->{ends}";
        say "lengths: file: $lengths; sequence: $location_info->{lengths}";
    }
    $n_tests++;
}

say "";
say "$tests_passed of $n_tests tests passed.";
if ($tests_passed == $n_tests) {
    say 'All tests passed';
} else {
    say 'Some tests failed. (' . ( $n_tests - $tests_passed ) . ' total)';
}
