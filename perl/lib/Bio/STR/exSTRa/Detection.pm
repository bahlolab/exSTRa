#!/usr/bin/env perl

=head1

within read STR finding

Some ideas to find STRs within reads. 

=cut

# Define a class to store our results
package Bio::STR::exSTRa::Detection;
use Data::Dumper;

use 5.014;
use Moose;
use MooseX::FollowPBP;
use namespace::autoclean;
use List::Util qw(max);

# Attributes
has 'repeat_unit' => (
    isa => 'Str',
    is  => 'rw',
    # trigger => sub {
    #     my $self = shift; 
    #     $self->clear_repeat_unit_canonical();
    # },
    trigger => \&_find_str,
    required => 1,
);

has 'sequence' => (
    isa => 'Str',
    is  => 'rw',
    trigger => \&_find_str,
    required => 1,
);

has 's_ru_rotation' => (
    isa => 'Str',
    is  => 'ro',
    writer => '_set_s_ru_rotation',
    default => '',
);

has 's_is_repeat_find' => (
    isa => 'Str',
    is  => 'ro',
    writer => '_set_s_is_repeat_find',
    default => '',
);

has 'is_repeat' => (
    isa => 'Str',
    is  => 'ro',
    writer => '_set_is_repeat',
    reader => 'get_is_repeat',
    default => '',
);

has 'min_bp_pre' => (
    isa => 'Int',
    is  => 'rw',
    default => 10,
    trigger => \&_find_str,
);

has 'max_bp_join' => (
    isa => 'Int',
    is  => 'rw',
    default => 18, # determined from known STRs, but may not be optimal or should really be dynamic (dependent on the lengths of found STRs)
    trigger => \&_find_str,
);

has 'min_bp' => (
    isa => 'Int',
    is  => 'rw',
    default => 15,
    trigger => \&_find_str,
);

has 'matches' => (
    isa => 'Bool',
    is  => 'ro',
    writer => '_set_matches',
);

has 'matches_same' => (
    isa => 'Str',
    is  => 'ro',
    writer => '_set_matches_same',
);

# Building
sub BUILD {
    my $self = shift;
    my $args = shift; 
    $self->{repeat_unit} = $args->{repeat_unit};
    $self->{sequence} = $args->{sequence};
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
    my $min_repeat_length_initial = $self->get_min_bp_pre; # in bp, adjust for best performance
    my $min_repeat_length = $self->get_min_bp; # in bp, adjust for best performance
    my $max_distance_length = $self->get_max_bp_join; # in bp, adjust for best performance

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
    # Now find the happy sequence
    # HMM?
    # Find sections that are pure repeats up to a threshold?
    # stitch long bits together?
    # threshold for mistakes? 1 repeat length? 2?
    my $i = 0;
    # initial repeat find
    my $min_repeat_starts_initial = $min_repeat_length_initial - $motif_len + 1;
    while($i < @motif_match) {
        # find first potential match
        if($motif_match[$i]) {
            # found first match
            my $start = $i;
            my $expected_next = $motif_match[$i] + 1;
            if($expected_next > $motif_len) {
                $expected_next = 1;
            }
            $i++;
            while(defined($motif_match[$i]) && $motif_match[$i] == $expected_next) {
                $expected_next++;
                if($expected_next > $motif_len) {
                    $expected_next = 1;
                }
                $i++;
            }
            if($i - $start + 1 >= $min_repeat_starts_initial) {
                # this appears to be the repeat
                foreach (@repeat_ok[($start..($i - 1))]) {
                    $_ = 1;
                }
            } else {
                # this appears to not be the repeat
                # do nothing
            }
            $i = max($start, $i - $motif_len); #make sure we don't miss the start of a repeat
        }
        $i++;
    }
    # join potential repeat regions
    my $repeat_ok_string = join ("", @repeat_ok);
    my $max_distance_starts_length = $max_distance_length + $motif_len;
    while($repeat_ok_string =~ /1(0{1,$max_distance_starts_length})1/gp) {
        my $centre = ${^MATCH};
        my $prem = ${^PREMATCH};
        my $postm = ${^POSTMATCH};
        while($centre =~ s/0/1/g) {}
        $repeat_ok_string = $prem . $centre . $postm;
    }

    # filter those that have repeats too short
    my $min_repeat_start_length = $min_repeat_length - $motif_len + 1;
    while($repeat_ok_string =~ /0(1{1,$min_repeat_start_length})0/gp) {
        # NOTE: repeats at the edge of the string will not be modified by this, but 
        # this is probably ok as we accept short repeats at the edges
        my $centre = ${^MATCH};
        my $prem = ${^PREMATCH};
        my $postm = ${^POSTMATCH};
        # Had to commnent this out, didn't appear to be sensible and caused an infinite loop
        #unless(length($prem) == 0 || length($postm) == 0) {
        #    # think twice about the threshold for the edges
        while($centre =~ s/1/0/g) {}
        #}
        $repeat_ok_string = $prem . $centre . $postm;
    }
    # say $sequence_raw;
    # say join ("", @motif_match);
    $self->_set_s_ru_rotation ( join ("", @motif_match) );
    # say join ("", @repeat_ok);
    # say $repeat_ok_string;
    $self->_set_s_is_repeat_find ( $repeat_ok_string );
    my $repeated_bases_ok = $repeat_ok_string;
    $repeated_bases_ok = $repeated_bases_ok . ("0" x ($motif_len - 1)); 
    # my $motif_len_m2 = $motif_len - 2;
    # my $repeat_ok_replace = (1 x $motif_len);
    # while($repeated_bases_ok =~ s/10.{$motif_len_m2}/$repeat_ok_replace/g) {}
    my @repeat_end_locations;
    while($repeated_bases_ok =~ /10/g) { # is the end is 1 of the $repeat_ok_string this will also make all of the end bases here 1
        push @repeat_end_locations, $-[0];
    }
    foreach (@repeat_end_locations) {
        substr $repeated_bases_ok, $_, $motif_len, (1 x $motif_len);
    }
    # say $repeated_bases_ok;
    $self->_set_is_repeat ( $repeated_bases_ok ); 
    return $repeated_bases_ok;
}

sub verify_repeat_found {
    # Check if the repeat location is found on the read
    # and set matches
    # Relies on the input sequence containing lower case letters 
    # for repeat bases and upper case for all others.
    my $self = shift;
    my $should_match = $self->get_sequence;
    while($should_match =~ s/[A-Z]/0/g) {}
    while($should_match =~ s/[a-z]/1/g) {}  
    my $result = $should_match eq $self->get_is_repeat || 0;
    $self->_set_matches ( $result );
    if($result) {
        $self->_set_matches_same ( 'equal' );
    } else {
        # give the bases that are unmatched
        my $bases_matching;
        my @shoulds = split //, $should_match;
        my @is_repeats = split //, $self->get_is_repeat;
        if (@shoulds == @is_repeats) {
            while(@is_repeats) {
                $bases_matching .= shift @is_repeats == shift @shoulds ? '1' : '0';
            }
            $self->_set_matches_same ( $bases_matching );
        } else {
            $self->_set_matches_same ( 'wrong_length' );
        }
    }
    return($result);
}

sub qualloc {
    # give a qualitative location for the STR
    # if given an argument then this should be a number corresponding to the number of bases
    my $self = shift;
    my $clipping = shift // 10;
    my $is_rep = $self->get_is_repeat;
    # warn "In qualloc, read clipping is ${clipping} bp\n";
    $is_rep = substr $is_rep, $clipping, length($is_rep) - 2 * $clipping; 
    for ($is_rep) {
        when (/^0+$/) { return 'outer' }
        when (/^0+1+0+$/) { return 'span' }
        when (/^1+0+$/) { return 'edge3p' }
        when (/^0+1+$/) { return 'edge5p' }
        when (/^1+$/) { return 'inner' }
        default { return 'NA' } 
    }
}

sub rep_in_read {
    # Give the bases of non-repeat, repeat then repeat in a read
    # Gives no result when the repeat locations cannot be decribed 
    # in this format such as multiple repeats, or when there is no
    # repeat in the read at all
    my $self = shift;
    my $clipping = shift // 10;
    my $is_rep = $self->get_is_repeat;
    if($is_rep =~ /^.{$clipping}(0*)(1+)(0*).{$clipping}$/) {
        return [length($1), length($2), length($3)];
    } else {
        return undef;
    }
}


sub location {
    # give the location of the STR with respect to the sequence
    my $self = shift;
    my $is_rep = $self->get_is_repeat;
    # determine if the repeat occurs more than once in the sequence
    my $here = 0;
    my @starts = ();
    my @ends = ();
    my @lengths = ();
    while($is_rep =~ s/^(0*)(1+)//) {
        my $len0 = length($1);
        my $len1 = length($2);
        push @starts, $here + $len0;
        push @lengths, $len1;
        $here += $len0 + $len1;
        push @ends, $here;
    }
    return { starts => \@starts, ends => \@ends, lengths => \@lengths };
}

no Moose; 

1;


