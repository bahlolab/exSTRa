#!/usr/bin/env perl

=head1

Within read STR scoring. Gives a score as to the number of times the repeat unit is seen within a read.

=cut

# Define a class to store our results
package Bio::STR::Score;
use Data::Dumper;

use 5.014;
use Moose;
use MooseX::FollowPBP;
use namespace::autoclean;
use List::Util qw(max sum);

# Types
# TODO: can I put these into a single file for both Bio::STR and Bio::STR::Score
use Moose::Util::TypeConstraints;
subtype 'Bio::STR::Score::Natural',
    as 'Int',
    where { $_ > 0 };
no Moose::Util::TypeConstraints;


# Attributes
has 'repeat_unit' => (
    isa => 'Str',
    is  => 'ro',
    # trigger => sub {
    #     my $self = shift; 
    #     $self->clear_repeat_unit_canonical();
    # },
    required => 1,
);

has 'sequence' => (
    isa => 'Str',
    is  => 'ro',
    required => 1,
);

has 's_ru_rotation' => (
    isa => 'Str',
    is  => 'ro',
    writer => '_set_s_ru_rotation',
    default => '',
);

has 'is_repeat' => (
    isa => 'Str',
    is  => 'ro',
    writer => '_set_is_repeat',
    reader => 'get_is_repeat',
    default => '',
);

has 'repeated_bases' => (
    isa => 'Int',
    is => 'ro',
    writer => '_set_repeated_bases',
);

has 'matchable_bases' => (
    # the maximum value that repeated_bases may take
    isa => 'Bio::STR::Score::Natural',
    is => 'ro',
    writer => '_set_matchable_bases',
);

# Building
sub BUILD {
    my $self = shift;
    $self->_find_str;
}

### Subroutines
sub _find_str {
    # finds the given STR in the string
    # Inputs: 
    #   motif => string of repeat motif
    #   sequence => query string of bases to find repeats
    # Optional inputs
    #   min_bp_pre => minimum number of required for finding potential repeat regions
    #   max_bp_join => the maximum distance away that two repeat regions can be if joined
    #   min_bp => minimum bp a repeat must be to be reported as a repeat
    my $self = shift;
    my $motif = $self->get_repeat_unit or die "No motif given";
    my $sequence_raw = $self->get_sequence or die "No sequence given";

    # The main work
    my $sequence = uc $sequence_raw;
    my $seq_len = length $sequence;
    my $motif_len = length $motif;
    unless ($motif =~ /^[ACGT]+$/) { die "Can only find exact DNA motifs" }
    unless ($sequence =~ /^[ACGTN]+$/) { die "Can only search within DNA" }
    my $motif_rota = $motif;
    my @motif_match = (0) x ($seq_len - $motif_len + 1);
    my @repeat_ok = @motif_match; # use for later
    for (my $rota_i = 1; $rota_i <= $motif_len; $rota_i++) {
        while($sequence =~ /(?=$motif_rota)/g) {
            $motif_match[$-[0]] = $rota_i;
        }
        $motif_rota =~ s/^(.)(.+)/$2$1/;
    }
    unless ($motif eq $motif_rota) { die "Program bug: should have fully rotated the repeat motif by now." } 
    $self->_set_s_ru_rotation ( join ("", @motif_match) );
    # Perform the scoring
    # Extend matches out
    my $extend = 0;
    push @motif_match, ((0) x ($motif_len - 1));
    foreach (@motif_match) {
        if($_ != 0) {
            # prepare to add more bases if need be
            $extend = $motif_len - 1;
        } elsif ($extend != 0) {
            $_ = 1;
            $extend--;
        }
    }
    my $repeated_bases = sum map { $_ > 0 } @motif_match;
    #$self->_set_matchable_bases ($seq_len - $motif_len + 1);
    $self->_set_matchable_bases ($seq_len);
    $self->_set_repeated_bases ( $repeated_bases ); 
    return;
}

sub repeated_proportion {
    # give the repeatd proportion
    my $self = shift;
    return ($self->repeated_bases / $self->matchable_bases);
}


no Moose; 

1;


