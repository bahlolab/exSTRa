#!/usr/bin/env perl

=head1

within read STR finding

Some ideas to find STRs within reads. 

=cut


use 5.014;
use warnings; 
use Data::Dumper;
use Bio::STR::exSTRa::Score;

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
    my $strscore = Bio::STR::exSTRa::Score->new (repeat_unit => $rep_unit, sequence => $seq);
    say $line[$headcols{Disease}];
    say $strscore->{sequence};
    say $strscore->repeated_bases . " / " . $strscore->matchable_bases;
    say "Score: " . $strscore->repeated_proportion;
    if($strscore->repeated_proportion <= 1) {
        $tests_passed++;
    } else {
        say "FAILED, score is >1";
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
