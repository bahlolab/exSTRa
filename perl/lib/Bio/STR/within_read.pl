#!/usr/bin/env perl

=head1

within read STR finding

Some ideas to find STRs within reads. 

=cut


use 5.014;
use warnings; 
use Data::Dumper;
use Bio::STR::Detection;

my %headcols;
my ($tests_passed, $n_tests) = (0, 0);
my $test_qualloc = '';
my $qualloc;
while (<>) {
    chomp;
    if(/^#END#/) {
        last; # don't bother with any more columns
    }
    my @line = split /\t/;
    if($. == 1) {
        @headcols{@line} = (0..$#line);
        if (defined ($headcols{qual_loc})) {
            $test_qualloc = 1;
        }
        next;
    }
    my ($rep_unit, $seq) = @line[@headcols{qw(repeat_unit sequence)}];
    if(!defined($rep_unit) || !defined($seq)) {
        die "There does not appear to be both of the columns 'repeat_unit' and 'sequence' in the input\n";
    }
    my $strdetect = Bio::STR::Detection->new (repeat_unit => $rep_unit, sequence => $seq);
    $strdetect->verify_repeat_found;
    if ($test_qualloc) {
        $qualloc = $strdetect->qualloc;
    }
    if($strdetect->matches && (!$test_qualloc || $qualloc eq $line[$headcols{qual_loc}])) {
        say $line[$headcols{Disease}] . ' matches';
        $tests_passed++;
    } else {
        say "\n" . $line[$headcols{Disease}] . ' does not match! Details:';
        for (qw(repeat_unit sequence matches_same)) {
            printf "%-12s %s\n", $_, $strdetect->{$_};
        }
        if ($test_qualloc) {
            printf "%-12s Expected = '%s', Detected = '%s'\n", 'qualloc', $line[$headcols{qual_loc}], $qualloc;
        }
        say "";

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
