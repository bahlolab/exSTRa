#
# svn $Revision: 1014 $
# svn $LastChangedDate: 2016-09-01 20:07:23 +1000 (Thu, 01 Sep 2016) $

=head1 Summary

Defines the STR object and its properties


=cut 

use 5.014;
use warnings;

package Bio::STR;

# use Spreadsheet::XLSX;
use autodie;
use Carp;

our(@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS, $VERSION);

use Exporter; 
$VERSION = 0.02;
@ISA = qw(Exporter); 

@EXPORT     = qw ();
@EXPORT_OK  = qw ();
%EXPORT_TAGS = ();

### Main code ###
{
package STR;

=head2
Bio::STR

Defines the properties of one simple tandem repeat

=cut

use Moose; 
use namespace::autoclean;

# Types
use Moose::Util::TypeConstraints;
subtype 'Natural',
    as 'Int',
    where { $_ > 0 };

subtype 'Natural0',
    as 'Int',
    where { $_ >= 0 };

subtype 'Strand',
    as 'Str',
    where { $_ =~ /^[-+]?$/ };

no Moose::Util::TypeConstraints;

# Attributes
has 'repeat_unit' => (
    isa => 'Str',
    is  => 'ro',
    trigger => sub {
        my $self = shift; 
        $self->clear_repeat_unit_canonical();
        $self->clear_repeat_unit_fwd();
    },
);

has 'repeat_unit_canonical' => (
    isa => 'Str',
    is  => 'ro', 
    lazy => 1, 
    clearer => 'clear_repeat_unit_canonical',
    default => sub { # To be inferred from repeat unit
        my $self = shift;  
        unless ($self->compound_strs_allowed) {
            $self->canonical_seq();
        }
    } 
);

has 'repeat_unit_fwd' => (
    # The forward repeat unit
    isa => 'Str',
    is  => 'ro',
    lazy => 1,
    clearer => 'clear_repeat_unit_fwd',
    default => sub {
        my $self = shift;
        if ($self->strand eq '-') {
            STR::_reverse_comp ($self->repeat_unit);
        } else {
            $self->repeat_unit;
        }
    }
);

# has 'location' => (
#     isa => '?',
#     is  => 'ro',
# );

has 'compound_strs_allowed' => (
    isa => 'Bool', 
    is => 'rw',
    default => 0,
);

has [qw(chrom name gene region)] => (
    isa => 'Str',
    is  => 'ro',
);

has [qw(strand)] => (
    isa => 'Strand',
    is  => 'ro', 
    default => '',
);

has [qw(start end)] => (
    isa => 'Natural',
    is  => 'ro',
    trigger => sub {
        # Check that start <= end
        my $self = shift;
        if (defined($self->start) && defined($self->end)) {
            $self->start <= $self->end or die "A STR was attempted to be loaded with start ($self->end) < end ($self->end) position. Chromosome (if defined): $self->chrom.";
        }
    }
);

has [qw(period consensus_size)] => (
    isa => 'Natural',
    is  => 'ro',
);

has [qw(per_match per_indel score entropy)] => (
    isa => 'Num',
    is  => 'ro',
);

has [qw(stable_repeats unstable_repeats)] => (
    isa => 'Str',
    is  => 'ro',
);

has 'read_loc_per_sample' => (
    isa => 'HashRef[STR::ReadLocation | STR::ReadLocation::ReadInspect]',
    is  => 'ro',
    default => sub { { } }, 
);

has 'read_loc_pooled' => (
    isa => 'HashRef',
    is  => 'ro',
    default => sub { { } }, 
);

has 'read_detect_size' => (
    isa => 'Int',
    is  => 'rw',
);

has 'read_detect_size_is_empirical' => (
    isa => 'Bool',
    is  => 'rw',
);

# Methods
sub start0 {
    # gives the 0-based coordinates
    my $self = shift;
    $self->start - 1;
}


sub canonical_seq {
    # takes a DNA seq, then gives a canonical form so that
    # it is the first alphabetically on both the forward and reverse strand
    my $self = shift;
    my $seq = uc $self->{repeat_unit};
    $seq =~ /^[ACGT]*$/ or warn "Sequence $seq does not appear to be DNA";
    my $seq_rev = scalar reverse $seq;
    $seq_rev =~ tr/ACGT/TGCA/;
    my $lowest = $seq;
    for(my $i=0; $i < length($seq); $i++) {
        my $foroard = (substr $seq, $i) . (substr $seq, 0, $i);
        my $reverse = (substr $seq_rev, $i) . (substr $seq_rev, 0, $i);
        if($foroard lt $lowest) { $lowest = $foroard };
        if($reverse lt $lowest) { $lowest = $reverse };
    }
    $lowest;
}


sub read_loc_calc {
    # Calculates read location of reads in the BAM files
    # for a single STR object
    # Usage: 
    # $str->read_loc_calc($bams);
    # $bams is a reference to a hash: keys=sample names, values=Bio::DB::Sam objects
    my $self = shift;
}

sub give_empirical_str_size {
    # Give the empirical STR size
    my $self = shift;
    my %inputs = @_;
    my $reference_fasta = $inputs{reference} or die "No FASTA file given";
    my $outer_distance = $inputs{outer_distance} // 800; # change depending on insert size
    # get sequence and repeat unit
    my $sequence = _get_sequence($reference_fasta, $self->{chrom}, $self->{start}, $self->{end});
    my $detect = Bio::STR::Detection->new (repeat_unit => $self->repeat_unit_fwd, sequence => $sequence); 
    my $locations = $detect->location();
    if(@{$locations->{lengths}} == 1) {
        $self->{read_detect_size} = $locations->{lengths}->[0];
    } elsif (@{$locations->{lengths}} == 0) {
        $self->{read_detect_size} = 0; # no location
    } else {
        $self->{read_detect_size} = ''; # multiple locations
    }
    $self->{read_detect_size_is_empirical} = 1;
    # TODO: test this
}

sub _reverse_comp {
    # gives the reverse complement of the input
    my $seq_rev = scalar reverse $_[0];
    $seq_rev =~ tr/ACGT/TGCA/;
    $seq_rev;
}

sub _get_sequence {
    # Check inputs
    my($reference, $chrom, $start, $end) = @_;
    $start =~ /^\d+$/ or die "Not given a positive integer ($start).";
    $end =~ /^\d+$/ or die "Not given a positive integer ($end).";
    my $sequence = qx/samtools faidx $reference $chrom:$start-$end/ or die "reading Fasta file failed";
    $sequence =~ s/^.+\n//; # Remove Fasta sequence description
    $sequence =~ s/\n//g; # Remove excess newlines
    return($sequence);
}

1; # Success
}

### Main code ###
{
package STR::DB;

=head2
Bio::STR::DB

Describes a collection of Bio::STR objects. Includes methods for importing these collections from files. 

=cut 

use 5.014;
use Carp;
use Moose; 
use namespace::autoclean;
use autodie;
use Spreadsheet::XLSX;
use Spreadsheet::ParseExcel;
use warnings;
#use STR;
use Bio::DB::Sam;
use Data::Dumper;
use Text::Iconv;
use Tie::IxHash;
use List::Util qw(all);

# Types
use Moose::Util::TypeConstraints;
subtype 'File',
    as 'Str',
    where { -e $_ },
    message { "This file ($_) does not appear to exist!" },
    ;

no Moose::Util::TypeConstraints;

# Attributes
has [qw(strs)] => (
    # Database of the STRs
    is  => 'ro',
    isa => 'HashRef[Bio::STR]',
    default => sub { {} },
);

has [qw(strs_by_name strs_by_gene)] => (
    # Database of the STRs
    is  => 'ro',
    isa => 'HashRef[Bio::STR]',
    default => sub { 
        my %hash; 
        tie %hash, 'Tie::IxHash';
        \%hash; 
    },
);

has 'fasta' => (
    is  => 'rw',
    isa => 'File',
    predicate => 'has_fasta',
);

has 'bams' => (
    # database of the bam files
    isa => 'HashRef[Bio::DB::Sam]',
    is  => 'ro',
    default => sub { {} },
);

has 'compound_strs_allowed' => (
    # Allows STRs that have more than one repeat if enabled 
    isa => 'Bool', 
    is => 'rw',
    default => 0,
);

has 'counts' => (
    isa => 'HashRef',
    is => 'ro',
    default => sub { {} },
);

# Methods
# Subs predeclare
sub _insert_to_hash(++);

sub read_str_database_UCSC_TRF {
    # Load a database of STRs from the Tandem Repeat Finder (Simple Repeats)
    # track of UCSC Genome Browser.
    my $self = shift; 
    my $str_data = shift;
    my $has_head = shift // 1; # second argument, true if the input file has a header, false if straight from UCSC without
    my $allow_multi = shift // ''; # third arg, multiple STR loci are allowed, except only the XXXXX is kept
    if(%{$self->strs} ne 0) {
        die "Trying to load STR database twice, not allowed.\n";
    }
    if ($str_data =~ /\.xlsx?$/) {
        die "This function does not read xls/xlsx files. Use read_str_database_excel instead. \n";
    }
    open STRDB, '<', $str_data;
    # Define columns by name
    my %coli;
    my %data_2_internal = (
        'chrom'         => 'chrom',
        'chromStart'    => 'start',
        'chromEnd'      => 'end',
        'period'        => 'period',
        'consensusSize' => 'consensus_size',
        'perMatch'      => 'per_match',
        'perIndel'      => 'per_indel',
        'score'         => 'score',
        'entropy'       => 'entropy',
        'sequence'      => 'repeat_unit', 
    );
    if($has_head) {
        my $head = <STRDB>;
        chomp $head;
        my @heads = split /\t/, $head;
        @coli{@heads} = (0..$#heads);
        # Verify required columns
        foreach (keys %data_2_internal) {
            exists($coli{$_}) or die "The column $_ is not in the STR input file $str_data.\n";
        }
    } else {
        @coli{qw(bin chrom chromStart chromEnd name period copyNum consensusSize 
            perMatch perIndel score A C G T entropy sequence)} 
            = (0..100);
    }
    # Import STRs
    while(<STRDB>) {
        chomp;
        my @line = split(/\t/);
        my $str = STR->new(compound_strs_allowed => $self->compound_strs_allowed);
        ++$line[$coli{'chromStart'}]; # account for 0-based coordinates
        foreach my $internal (keys %data_2_internal) {
            $str->{$data_2_internal{$internal}} = $line[$coli{$internal}];
        }
        $str->repeat_unit_canonical; # Generate this value as it's normally lazy
        my $key = "$str->{chrom}:$str->{start}-$str->{end}:$str->{repeat_unit}";
        if(defined($self->strs->{$key})) {
            unless($allow_multi) {
                die "Repeated location of repeat at location $key";
            }
            next; 
        }
        $str->{name} = $key; # set name for STR
        $self->strs->{$key} = $str; # Add STR to hash
        
    }
    close STRDB;
    my $strs_imported = scalar(keys %{$self->strs});
    warn "STRs imported: $strs_imported\n";
    unless($strs_imported) {
       die "No STRs were imported";
    }
}

sub read_str_database_excel {
    # Load a database of STRs, we should be able to accept multiple formats.
    # This method will accept multiple tab delimted inputs, provided the column
    # names are specified.
    # Usage: 
    #   $strdb->read_str_database_tab("file.xlsx", %feature_names);
    # where %feature names has key=attribute name, value=header name. start0 may also be used
    my $self = shift; 
    my $excel_file = shift;
    my %feature_names = @_;
    my $read_detect_in_file = defined($feature_names{read_detect_size});
    if(%{$self->strs} ne 0) {
        die "Trying to load STR database twice, not allowed.\n";
    }

    my $excel;
    if ($excel_file =~ /\.xlsx$/) {
        my $converter = Text::Iconv -> new ("utf-8", "windows-1251");
        $excel = Spreadsheet::XLSX -> new ($excel_file, $converter); 
    } elsif ($excel_file =~ /\.xls$/) {
        my $parser   = Spreadsheet::ParseExcel->new();
        $excel = $parser->parse($excel_file);
    } else {
        die "Cannot determine spreadsheet format by extension (expected to see .xls or .xlsx). File $excel_file\n";
    }

    foreach my $worksheet (@{$excel -> {Worksheet}}) {
        my ( $row_min, $row_max ) = $worksheet->row_range();
        my ( $col_min, $col_max ) = $worksheet->col_range();
        # read header
        my @heads;
        foreach my $col (0 .. $col_max) {
                my $cell = $worksheet->get_cell( $row_min, $col );
                push @heads, ($cell ? $cell->value() : '');
        }
        my %cols;
        @cols{@heads} = (0..$#heads);
        unless (all {defined} @cols{(values %feature_names)}) {
            warn "When importing STRs for file $excel_file,\nsheet $worksheet->{Name} skipped due to not containing all specified columns.\n";
            warn "This may be ok if there are multiple sheets in the Excel file.\n";
            warn "Input column names and column number:\n";
            warn ((Dumper \%cols). "\n");
            warn "Bio::STR class attribute names and names we expect to find in Excel file:\n";
            warn ((Dumper \%feature_names). "\n");
            next;
        }
        # Read body
        for my $row ( ($row_min + 1) .. $row_max ) {
            my $str = STR->new(compound_strs_allowed => $self->compound_strs_allowed);
            for my $feat (keys %feature_names) {
                my $cell = $worksheet->get_cell( $row, $cols{$feature_names{$feat}} );
                my $cvalue;
                if ($cell) {
                    $cvalue = $cell->value();
                    if($feat eq 'start0' && $cvalue =~ /^\d+$/) {
                        $str->{start} = $cvalue + 1;
                    } else {
                        $str->{$feat} = $cvalue;
                    }
                } 
            }
            # Check STR has a location
            unless (defined($str->{chrom}) && defined($str->{start}) && defined($str->{end}) && defined($str->{repeat_unit})) {
                warn "When importing STRs for $excel_file,\nsheet $worksheet->{Name}, row " . ($row + 1) . " skipped.\n";
                no warnings; 
                warn "Values found that were required were: Chrom: $str->{chrom}, start: $str->{start}, end: $str->{end}, repeat_unit: $str->{repeat_unit}\n";
                use warnings; 
                next;
            }
            # mark read detect info as emperical if required
            if($read_detect_in_file) {
                $str->{read_detect_size_is_empirical} = '';
            }
            # Write the STR to the hash
            my $key = "$str->{chrom}:$str->{start}-$str->{end}:$str->{repeat_unit}";
            die "Repeated location of repeat at location $key" if defined($self->strs->{$key});
            $str->repeat_unit_canonical; # Generate this value as it's normally lazy
            $self->strs->{$key} = $str; # Add STR to hash
            if (defined($str->{name})) { $self->strs_by_name->{$str->{name}} = $str; }
            if (defined($str->{gene})) { $self->strs_by_gene->{$str->{gene}} = $str; }
        }
     }
     my $strs_imported = scalar(keys %{$self->strs});
     warn "STRs imported: $strs_imported\n";
     unless($strs_imported) {
        die "No STRs were imported";
     }
}

sub keep_repeat_unit_sizes {
    # Filters repeat database hash by the repeat unit size
    # Usage: 
    # $strdb->filter_repeat_unit_size(MIN, MAX)
    my $self = shift; 
    my $min = shift // 2;
    my $max = shift // 6; 
    if ($min > $max || $min !~ /^\d+$/ || $max !~ /^\d+$/) {
        die "Input to keep_repeat_unit_sizes invalid, require 2 positive integers where MIN <= MAX. MIN=$min, MAX=$max"
    }
    warn "Only keeping repeats with units $min to $max bp.\n";
    foreach my $key (keys %{$self->strs}) {
        my $cs = $self->strs->{$key}->{consensus_size};
        if ($cs < $min || $cs > $max) {
            delete $self->strs->{$key};
        }
    }
    warn "Now have " . scalar(keys %{$self->strs}) . " remaining STRs\n";
}

sub read_bams {
    # Reads in BAM files to be processed
    # Usage: 
    # $strdb->read_bams(\%bam_filenames, @options)
    # where %bams has keys=sample names, values = bam file paths
    my $self = shift;
    my $bam_filenames = shift;
    unless($self->has_fasta) { croak "Require FASTA index to be specified before reading in BAMs" }
    foreach my $sample (keys %{$bam_filenames}) {
        warn "Loading BAM for sample $sample.\n";
        if(exists($self->bams->{$sample})) { warn "Reloading existing loaded BAM file for sample $sample.\n" }
        my $bam = Bio::DB::Sam->new(
            -bam => $bam_filenames->{$sample},
            -fasta => $self->fasta,
            -autoindex => 1,
            -expand_flags => 1,
            @_,
        );
        $self->bams->{$sample} = $bam;
    }
}

sub read_bams_array {
    # Reads in BAM files to be processed from an array, deducing 
    # sample names from the BAM files themselves
    # Usage: 
    # $strdb->read_bams(\@bam_filenames, @options)
    # where %bams has keys=sample names, values = bam file paths
    my $self = shift;
    my $bam_filenames = shift;
    for my $bam_file (@$bam_filenames) {
        my $bam = Bio::DB::Sam->new(
            -bam => $bam_file,
            -autoindex => 1,
        );
        $bam->header->text =~ /\tSM:([^\t]+)\t/;
        my $sample_name = $1;
        $self->read_bams({ $sample_name => $bam_file }, @_);
    }
}

sub assess_str_reads {
    my $self = shift;
    $self->assess_str_reads_by_readinspect (@_);
}

sub assess_str_reads_by_alignment {
    # Gives location of reads for each sample for each STR by reading the alignments
    # Usage
    # $strdb->assess_str_reads_by_alignment(DIST);
    # where DIST is a number of the outer distance to look for reads
    my $self = shift;
    # TODO: give this a named hash input
    my $outer_distance = shift // 100;
    if(%{$self->bams} eq 0) {
        die "No BAMs have been loaded yet in assess_str_reads.";
    }
    foreach my $str (values %{$self->strs}) {
        foreach my $sample (keys %{$self->bams}) {
            $str->read_loc_per_sample->{$sample} = STR::ReadLocation->new( 
                bam => $self->bams->{$sample}, 
                str => $str, 
                outer_distance => $outer_distance, 
            );
        }
    }
}

sub assess_str_reads_by_readinspect {
    # Gives location of reads for each sample for each STR with the readinspect algorithm
    # Usage
    # $strdb->assess_str_reads_by_readinspect(DIST);
    # where DIST is a number of the outer distance to look for reads
    my $self = shift;
    my %inputs = @_;
    my $outer_distance = $inputs{outer_distance} // 800; # should change depending on the average insert sizes
    my $read_trim_static = $inputs{read_trim_static} // 10; # how many bp to trim reads from both ends in qualloc
    my $print_only = $inputs{print_only} // 0;
    my $print_read_name = $inputs{print_read_name} // 1;
    warn "The options to assess_str_reads_by_readinspect are:\n";
    warn ((Dumper \%inputs) . "\n");
    if(%{$self->bams} eq 0) {
        die "No BAMs have been loaded yet in assess_str_reads.";
    }
    if(1) {
        warn "Read trim static is ${read_trim_static}bp\n";
    }
    if ($print_only) {
        if(defined($inputs{give_rep_in_read}) && $inputs{give_rep_in_read}) {
            print join ("\t", qw(locus sample a b c strand));
        } elsif (defined($inputs{give_score}) && $inputs{give_score}) {
            print join ("\t", qw(locus sample rep mlength));
        }
        if($print_read_name) {
            print "\tread_id";
        }
        say '';
    }
    foreach my $str (values %{$self->strs}) {
        foreach my $sample (keys %{$self->bams}) {
            # here instead just directly print out each output
            my $strloc_ri = STR::ReadLocation::ReadInspect->new( 
                bam => $self->bams->{$sample}, 
                str => $str, 
                outer_distance => $outer_distance, 
                read_trim_static => $read_trim_static,
                %inputs, # pass other inputs to create new object
            );
            if($print_only) {
                if(defined($inputs{give_rep_in_read}) && $inputs{give_rep_in_read}) {
                    foreach my $rir (@{$strloc_ri->rep_in_read}) {
                        say join ("\t", (
                            _disease_shorten($str->name),
                            $sample,
                            $rir->a,
                            $rir->b,
                            $rir->c,
                            $rir->strand,
                            ($print_read_name ? $rir->bam_read->name : ()),
                        ));
                    }
                } elsif (defined($inputs{give_score}) && $inputs{give_score}) {
                    foreach my $scr (@{$strloc_ri->score}) {
                        say join ("\t", (
                            _disease_shorten($str->name),
                            $sample,
                            $scr->repeat_starts,
                            $scr->mlength,
                            ($print_read_name ? $scr->bam_read->name : ()),
                        ));
                    }
                } else {
                    croak "Instant printing of this kind not yet implemented";
                }
            } else {
                $str->read_loc_per_sample->{$sample} = $strloc_ri;
            }
        }
    }
}

sub assess_reference_readinspect {
    # Gives the size of STRs in bp as determined by the readinspect algorithm
    # Also gives potential warnings for multiple occurances of the repeat motif
    # for each locus for the given DIST size.
    # Usage
    # $strdb->assess_reference_readinspect(DIST)
    # DIST: maximum distance to search for similar repeats
    my $self = shift;
    my $outer_distance = shift // 800; # should change depending on the average insert sizes
    unless ($self->has_fasta) {
        die "STR DB requires a FASTA file to be specified";
    }
    # STR has read_detect_size
    foreach my $str (values %{$self->strs}) {
        $str->give_empirical_str_size(reference => $self->{fasta}, outer_distance => $outer_distance)
        
    }
}


sub give_insert_sizes {
    # Gives insert sizes of the BAMs
    # Usage: 
    # give_insert_sizes(seqID, start, end)
    my $self = shift; 
    my $chrom = shift // 'chr2';
    my $start = shift // 1;
    my $end = shift // 10000000;
    warn "GIS";
    foreach my $key (keys %{$self->bams}) {
        warn "New BAM";
        my $bam = $self->bams->{$key};
        warn "BAM identified";
        my @pairs = $bam->features(
            -type => 'read_pair', 
            -seq_id => $chrom,
            -start => $start, 
            -end   => $end,
        );
        warn "Alignments loaded";
        foreach my $pair (@pairs) {
            # warn Dumper($pair);
            my $length                    = $pair->length;   # insert length
            my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
            #warn Dumper($first_mate);
            #warn Dumper($second_mate);
            my $f_start = $first_mate->start;
            my $s_start = defined($second_mate) ? $second_mate->start : "NA" ;
            my $f_end = $first_mate->end;
            my $s_end = defined($second_mate) ? $second_mate->end : "NA" ;
            say "$length\t$f_start\t$f_end\t$s_start\t$s_end";
        }
    }

}

sub str1 {
    # Gives an arbitrary STR to access general STR properties
    my $self = shift; 
    my (undef, $value) = each %{$self->strs};
    # keys %{$self->strs} # No need to reset iterator I think
    $value;
}

sub all_counts_give {
    # generates a table of raw counts
    my $self = shift; 
    # TODO: edit here
    my %all_counts;
    tie(%all_counts, 'Tie::IxHash');
    for my $sample (sort keys %{$self->bams}) {
        warn "Analysing $sample\n";
        my $bam = $self->bams->{$sample};
        my @haplos = ('reference1', ($sample =~ /expanded/ ? 'expanded' : 'reference2'));
        my $all_counts_1 = _give_all_counts($self, $sample, \@haplos);
        $all_counts{$sample} = $all_counts_1;
    }
    $self->{counts} = \%all_counts;
}

sub _give_all_counts {
    # Giving all the different kinds of counts, in a hash
    my ($strs, $sample, $haplos) = @_;
    my %all_counts_1;
    tie(%all_counts_1, 'Tie::IxHash');
    my $naming_scheme = %{$strs->strs_by_name} ? 'strs_by_name' : 'strs';
    for my $str_name (keys %{$strs->{$naming_scheme}}){
        my %counts_1_str_name;
        tie(%counts_1_str_name, 'Tie::IxHash');
        my $c_ref;
        my %count_h;
        tie(%count_h, 'Tie::IxHash'); # @count_h is ordered
        # spanning counts
        for my $count_type (qw(up pe)) {
            my $method = "counts_$count_type";
            $c_ref = $strs->{$naming_scheme}->{$str_name}->read_loc_per_sample->{$sample}->$method();
            %count_h = ();
            @count_h{map { "${count_type}_$_" } (keys %$c_ref)} = values %$c_ref;
            _insert_to_hash(%counts_1_str_name, %count_h);
        }
        # Add to the hash
        $all_counts_1{$str_name} = \%counts_1_str_name;
    }
    # Return value
    \%all_counts_1;
}

sub _insert_to_hash(++) {
    my $hash = shift;
    my $insert = shift; 
    while(my ($key, $value) = each %$insert) {
        $hash->{$key} = $value;
    }
}

sub all_counts_print {
    # Print all those counts into a tab delimted table
    my $self = shift;
    my $all_counts = $self->counts;
    my $str_name_col_head = 'locus'; #TODO: should this depend on whether this is a disease or not?
    # Print header
    my @datacols = keys ((each ((each $all_counts)[1]))[1]);
    # Header: sample, locus, data 
    say join("\t", ('sample', $str_name_col_head, @datacols));
    
    # Print body
    for my $sample (keys %$all_counts) {
        my $sample_counts = $all_counts->{$sample};
        for my $disease (keys %$sample_counts) {
            my $disease_counts = $sample_counts->{$disease};
            say join("\t", ($sample, _disease_shorten($disease), @{$disease_counts}{ @datacols })); 
        }
    }
    return 1;
}

sub rep_in_read_print {
    # Print the repeats in read counts to a tab delimited file
    my $self = shift;
    my $str_name_col_head = 'locus'; #TODO: should this depend on whether this is a disease or not?
    my %inputs = @_;
    my $read_name = $inputs{read_name} // 0;
    # Format could be:
    # locus, sample, a, b, c, strand? (once per read)
    # optional read ID column?
    say join("\t", (qw(locus sample a b c strand), ($read_name ? 'read_id': ())));
    foreach my $str (values %{$self->strs}) {
        foreach my $sample (keys %{$self->bams}) {
            foreach my $rir (@{$str->read_loc_per_sample->{$sample}->rep_in_read}) {
                say join ("\t", (
                    _disease_shorten($str->name),
                    $sample,
                    $rir->a,
                    $rir->b,
                    $rir->c,
                    $rir->strand,
                    ($read_name ? $rir->bam_read->name : ()),
                ));
            }
        }
    }
}

sub _disease_shorten {
    # shorten a disease name
    map { if(/\((\w+)\)/) { $1 } else { $_ } } @_;
}

1; # Success
}

### Main code ###
{
package STR::ReadLocation;
use Carp;
use Moose; 
use namespace::autoclean;
use autodie;
#use STR;
use Bio::DB::Sam;
use Data::Dumper;

use constant REGIONSSE => qw (
    00
    01
    11
    02
    12
    22
);
    
use constant REGIONSPE => qw (
    0000
    0001
    0002
    0011
    0012
    0022
    0100
    0101
    0102
    0111
    0112
    0122
    0201
    0202
    0212
    0222
    1101
    1111
    1112
    1122
    1211
    1212
    1222
    2212
    2222
    00xx
    01xx
    11xx
    02xx
    12xx
    22xx
);

# ReadLocation->new( bam => $bam_object, str => $str_object );

#
# Types


# Object contruction
sub BUILD {
    # TODO add read location determining methods
    my $self = shift;
    $self->determine_location_se();
    $self->determine_location_pe();
}


# Attributes
has 'bam' => (
    # given bam file
    isa => 'Bio::DB::Sam',
    is  => 'ro', 
    required => 1,
);

has [qw(outer_distance)] => (
    isa => 'Natural', 
    is  => 'ro',
    default => 100,
);

has 'str' => (
    isa => 'STR',
    is  => 'ro',
    required => 1,
);

has [REGIONSSE] => (
    # Single end reads, each group is mutally exclusive
    is  => 'ro',
    isa => 'ArrayRef[Bio::DB::Bam::Alignment]',
    default => sub { [] },
);

has [REGIONSPE] => (
    # Paired end reads, each group is mutally exclusive
    is  => 'ro',
    isa => 'ArrayRef[Bio::DB::Bam::Alignment]',
    default => sub { [] },
);

has 'regions_se_names' => (
    isa => 'ArrayRef[Str]',
    is => 'ro', 
    lazy => 1,
    default => sub { [REGIONSSE] },
);

has 'regions_pe_names' => (
    isa => 'ArrayRef[Str]',
    is => 'ro', 
    lazy => 1,
    default => sub { [REGIONSPE] },
);

has 'rep_in_read' => (
    isa => 'ArrayRef[STR::Rep_in_read]',
    is => 'ro',
    default => sub { [] },
);

has 'score' => (
    isa => 'ArrayRef[STR::Score]',
    is => 'ro',
    default => sub { [] },
);


# Methods
sub determine_location_se {
    # Categorise single ends of reads into their location relative to STRs
    my $self = shift;
    my @reads = $self->bam->get_features_by_location(-seq_id => $self->str->{chrom},
                                                 -start  => $self->str->{start} - $self->outer_distance,
                                                 -end    => $self->str->{end} + $self->outer_distance
                                                 );
    my $str_start = $self->str->start;
    my $str_end = $self->str->end;
    READPAIRS: for my $r (@reads) {
        # Categorise each read into where it lies over the STR
        if(_quality_filter($r, undef, 3)) {
            next READPAIRS;
        }
        my $start = $r->start;
        my $end = $r->end;
        my $key = '';
        for my $pos ($start, $end) {
            if ($pos < $str_start) {
                $key .= '0';
            } elsif ($pos <= $str_end) {
                $key .= '1';
            } else {
                $key .= '2';
            }
        }
        if (!exists($self->{$key})) { warn "No paired value for key $key." };
        push @{$self->{$key}}, $r;
    }
}

sub counts_se {
    # Give counts of the alignments with respect to single-ends
    # TODO: This should be able to use s_all???
    my $self = shift; 
    my %counts;
    tie(%counts, 'Tie::IxHash');
    for my $region (REGIONSSE) {
        $counts{$region} = scalar(@{$self->$region});
    }
    \%counts;
}

sub counts_up {
    # Synonym for counts_se
    counts_se(@_);
}

sub determine_location_pe {
    # Categorise paired end reads into their location relative to STRs
    my $self = shift;
    my @pairs = $self->bam->get_features_by_location(
                                                 -type   => 'read_pair',
                                                 -seq_id => $self->str->{chrom},
                                                 -start  => $self->str->{start} - $self->outer_distance,
                                                 -end    => $self->str->{end} + $self->outer_distance
                                                 );
    my $str_start = $self->str->start;
    my $str_end = $self->str->end;
    READPAIRS: for my $pair (@pairs) {
        # Categorise each paried-read into its overlap with the STR
        my ($first_mate, $second_mate) = $pair->get_SeqFeatures;
        if(_quality_filter($first_mate, $second_mate, 3)) {
            next READPAIRS;
        }
        my $start1 = $first_mate->start;
        my $end1 = $first_mate->end;
        my ($start2, $end2);
        # TODO: INSERT filters for reads, maybe as a function so we don't copy paste
        if(!defined($second_mate)) { 
            $start2 = 'x';
            $end2 = 'x';
            unless($start1 < $end1) {
                 warn "Read location is not as expected and being skipped. \$start1=$start1, \$end1=$end1, unpaired read.";
                 next;
            }
        } else { 
            $start2 = $second_mate->start;
            $end2 =  $second_mate->end;
            unless($start1 < $end1 && $start1 <= $start2 && $start2 < $end2 && $end1 <= $end2) {
                 warn "Read location is not as expected and may cause errors. \$start1=$start1, \$end1=$end1, \$start2=$start2, \$end2=$end2.";
                 # TODO: Deal with these reads better
            }
        }
        my $key = '';
        for my $pos ($start1, $end1, $start2, $end2) {
            if ($pos eq 'x') {
                $key .= 'x';
            } elsif ($pos < $str_start) {
                $key .= '0';
            } elsif ($pos <= $str_end) {
                $key .= '1';
            } else {
                $key .= '2';
            }
        }
        if (!exists($self->{$key})) { warn "No paired value for key $key." };
        push @{$self->{$key}}, $pair;
    }
}

sub counts_pe {
    # Give counts of the alignments with respect to paired-ends
    my $self = shift; 
    my %counts;
    tie(%counts, 'Tie::IxHash');
    for my $region (REGIONSPE) {
        $counts{$region} = scalar(@{$self->$region});
    }
    \%counts;
}


sub s_all {
    # Gives all the single reads around and in the STR
    my $self = shift;
    map { @{$self->{$_}} } REGIONSSE;
}

sub s_all_out {
    # Gives all the single end reads completely outside of the STR,
    my $self = shift;
    (@{$self->{00}}, @{$self->{11}});
}

sub s_all_over {
    # Gives all the single end reads with any overlap of the STR (including those inside)
    my $self = shift;
    (@{$self->{01}}, @{$self->{02}}, @{$self->{11}}, @{$self->{12}});
}

sub p_all {
    # Gives all the paired-end reads around and in the STR
    my $self = shift;
    map { @{$self->{$_}} } REGIONSPE;
}

sub p_all_over {
    # Gives all the paired-end reads that are in any way inside the STR, 
    # including those that are fully within the STR.
    my $self = shift;
    my @strs_out; 
    for my $code (REGIONSPE) {
        $code =~ /^(\d)\d*(\d)x*$/ or die "BUG";
        if ($1 == 1 || $2 == 1 || ($1 == 0 && $2 == 2)) {
            # Add these reads to those that are over the STR
            if (!defined($self->{$self->{$code}})) {
                die "INTERNAL BUG 697: Issues with $code, have no attribute " . $self->{$code};
            } else {
                push @strs_out, @{$self->{$self->{$code}}};
            }
        }
    }
    @strs_out;
}

sub p_spanning {
    # Gives all the paired-end reads that span the STR. Note this is for paired-end reads only
    my $self = shift;
    my @strs_out; 
    for my $code (REGIONSPE) {
        $code =~ /^(\d)\d*(\d)x*$/ or die "BUG";
        if ($1 == 0 && $2 == 2) {
            # Add these reads to those that are over the STR
            if (!defined($self->{$self->{$code}})) {
                die "INTERNAL BUG 727: Issues with $code, have no attribute " . $self->{$code};
            } else {
                push @strs_out, @{$self->{$self->{$code}}};
            }
        }
    }
    @strs_out;
}

sub _quality_filter {
    # provides a duplicate and mapping quality filter for pairs of reads, or a single end if the second read is undefined
    my ($first_mate, $second_mate, $qual_threshold) = @_;
    if($first_mate->get_tag_values('DUPLICATE')) {
        if(defined($second_mate)) {
            if($second_mate->get_tag_values('DUPLICATE')) {
                # Both reads marked as a duplicate
                return 1;
            } else {
                warn "WARNING: only first (chromosome position) read is marked as duplicate for read ". $first_mate->name ."\n";
            }
        } else {
            # The one mapped read has been marked as a duplicate
            return 1;
        }
    } elsif (defined($second_mate) && $second_mate->get_tag_values('DUPLICATE')) {
        warn "WARNING: only second (chromosome position) read is marked as duplicate for read ". $second_mate->name ."\n";
    }
    # Filter out reads with low mapping quality for both ends:
    #TODO: make $qual_threshold be setable from outside module
    if($first_mate->qual < $qual_threshold && (!defined($second_mate) || $second_mate->qual < $qual_threshold)) {
        # Both reads have a low mapping quality, so filter as they are unreliable
        return 1;
    }
    return 0;
}

}

### class for the location of reads using read inspection ###
{
package STR::ReadLocation::ReadInspect;
# A package to replace the STR::ReadLocation class, instead deriving counts
# from the Bio::STR::detection module/class to determine STR location with
# the read content. 
use Carp;
use Moose;
use namespace::autoclean;
extends 'STR::ReadLocation';
use Bio::STR::Detection;
use Bio::STR::Score;
# use Bio::DB::Bam::AlignWrapper;

# Attributes of STR::ReadLocation::ReadInspect

has [qw(read_trim_static)] => (
    isa => 'Natural0', 
    is  => 'ro',
    default => 10,
);

has [qw(give_qualloc give_rep_in_read give_score)] => (
    isa => 'Bool',
    is  => 'ro',
    default => 0,
);

# Methods of STR::ReadLocation::ReadInspect

sub determine_location_se {
    # Categorise single ends of reads into their location relative to STRs
    my $self = shift;
    my %inputs = @_;
    my $read_trim_static = $self->read_trim_static; # how many bp to trim reads from both ends in qualloc
    my $debug = '';
    # Check that this will actually do something
    unless ($self->give_qualloc || $self->give_rep_in_read || $self->give_score) {
        die "determine_location_se for read inspect requires either give_qualloc, give_rep_in_read or give_score to be turned on";
    }
    if($self->give_score && $read_trim_static != 0) {
        warn "Warning: Trimming is not implemented for give_score.\n"; 
    }
    if (defined($self->str->{name})) {
        warn 'Working on the locus ' . $self->str->{name} . "\n"; # TOREMOVE / on verbose
    } else {
        warn 'Working on the locus at ' . $self->str->{chrom} . ":" . $self->str->{start} . 
            "-" . $self->str->{end} . "\n"; # TOREMOVE / on verbose
    }
    #$debug = 1; # show extra info for some reads
    my %target_reads;
    my @target_reads_array = ();
    # @target_reads_array = qw (expanded_SCA10-1029 expanded_SCA10-1057 expanded_SCA10-21 expanded_SCA10-215 expanded_SCA10-233 expanded_SCA10-241 expanded_SCA10-271 expanded_SCA10-283 expanded_SCA10-301 expanded_SCA10-345 expanded_SCA10-349 expanded_SCA10-429 expanded_SCA10-449 expanded_SCA10-455 expanded_SCA10-471 expanded_SCA10-483 expanded_SCA10-491 expanded_SCA10-545 expanded_SCA10-55 expanded_SCA10-571 expanded_SCA10-637 expanded_SCA10-645 expanded_SCA10-665 expanded_SCA10-67 expanded_SCA10-679 expanded_SCA10-857 expanded_SCA10-929 expanded_SCA10-941 expanded_SCA10-999 );
    #@target_reads_array = qw (reference2_SCA3-511);
    #@target_reads{@target_reads_array} = (1) x @target_reads_array;
    # NOTE: simulations have no unmapped reads, need to check we retrieve these correctly
    my $pad_distance = 1000; # To try to ensure we get proper pairs
    my $min_start = $self->str->{start} - $self->outer_distance; # TODO: we need to make sure we skip over sites with no location
    my $max_end = $self->str->{end} + $self->outer_distance;
    my @pairs = $self->bam->get_features_by_location(
                                                 -type   => 'read_pair',
                                                 #-type   => ['read_pair', 'match'],
                                                 -seq_id => $self->str->{chrom},
                                                 -start  => $min_start - $pad_distance,
                                                 -end    => $max_end + $pad_distance,
                                                 );
    my $str_start = $self->str->start;
    my $str_end = $self->str->end;
    my $rep_unit = $self->str->repeat_unit_fwd;
    my %STRdetection2code_forw_strand = (
        'span' => '02',
        'inner' => '11',
        'edge5p' => '01',
        'edge3p' => '12',
        'outer' => '22', # we can't actually ever properly know this with these methods, could be dangerous to actually use
    );
    # my %STRdetection2code_rev_strand = (
    #     %STRdetection2code_forw_strand, 
    #     'edge5p' => '12',
    #     'edge3p' => '01',
    #     'outer' => '00', # we can't actually ever properly know this with these methods, could be dangerous to actually use
    # );
    READPAIRS: for my $pair (@pairs) {
        # Categorise each paried-read into its overlap with the STR
        # Directions are on forward strand
        my ($first_mate, $second_mate) = $pair->get_SeqFeatures;
        # TODO: Check that we don't have an unmapped read pair that we will be missing (make these statements earlier) [Is this still to do Rick?]
        # TODO: test these two filters: 
        # Filter read marked as duplicate:
        if(STR::ReadLocation::_quality_filter($first_mate, $second_mate, 3)) {
            next READPAIRS;
        }
        ANCHORBYMATE: for my $i (0, 1) {
            # set the first and second mates
            my $anchor_mate = $i ? $second_mate : $first_mate;
            my $detect_mate = $i ? $first_mate : $second_mate;
            if(!defined($anchor_mate) || $anchor_mate->get_tag_values('UNMAPPED')) {
                # We don't want to look at anchors that we don't have or are unmapped
                next ANCHORBYMATE;
            }
            my $interesting_read = '';
            if ($debug && $target_reads{$anchor_mate->name}) {
                warn "\nThis is an interesting read, " . $anchor_mate->name . "\n";
                $interesting_read = 1;
            }
            my $start1 = $anchor_mate->start;
            my $end1 = $anchor_mate->end;
            # check anchor read,
            my $anchor_direction = 0;
            if($start1 >= $min_start && $end1 < $str_end && ! $anchor_mate->get_tag_values('REVERSED')) {
                # The direction of the read is towards 3' ( ----> ).
                $anchor_direction = 1; # 3'
            } elsif ($end1 <= $max_end && $start1 > $str_start && $anchor_mate->get_tag_values('REVERSED')) {
                # The direction of the read is towards 5' ( <---- ).
                $anchor_direction = -1; # 5'
            }
            if ($anchor_direction) {
                # The start is as we would expect for an STR on the mate.
                # We have been very tolerant in this implementation.
                my $mate_recovered = 0;
                if(!defined($detect_mate)) {
                    # we want to ensure we have the mate of the read, it wasn't loaded so attempt to
                    # find it in the BAM file
                    $detect_mate = _find_mate ($anchor_mate);
                    unless(defined($detect_mate)) {
                        next ANCHORBYMATE;
                    }
                    $mate_recovered = 1;
                }
                # Check the sequence has been reversed if required
                my $sequence = $detect_mate->query->dna;
                if ($anchor_direction == 1 xor $detect_mate->get_tag_values('REVERSED')) {
                    # correct for when this read has been unexpectedly not reversed, maybe because it was un/mismapped
                    $sequence = STR::_reverse_comp $sequence;
                    if ($interesting_read) {
                        warn 'The mate was unexpectedly reversed ' . $detect_mate->name . "\n";
                    }
                }
                # Filter mates that are mapped such that they do not overlap the STR (and doesn't appear to be due to mismapping)
                # TODO: filter reads whose mapping suggests they are entirely not within the repeat.
                # TODO: Allow option to filter more harshly
                # TODO: think about +-1 differences (how have I defined end)
                #  $detect_mate->
                # TODO: don't run filter on recovered reads
                my $seq_len = length($sequence);
                unless ($mate_recovered || 
                    ($detect_mate->start + $seq_len > $str_start && $detect_mate->end - $seq_len < $str_end + 1)
                    ) {
                    next ANCHORBYMATE;
                }
                # Score
                if ($self->give_score) {
                    my $detect = Bio::STR::Score->new (repeat_unit => $rep_unit, sequence => $sequence); 
                    my $score = STR::Score->new (
                                repeat_starts => $detect->repeated_bases,
                                mlength => $detect->matchable_bases,
                                bam_read => $detect_mate,
                            );
                    if (defined($score)) {
                        push @{$self->score}, $score;
                    } 
                } 
                # Qual loc or rep in read
                if ($self->give_qualloc || $self->give_rep_in_read) {
                    my $detect = Bio::STR::Detection->new (repeat_unit => $rep_unit, sequence => $sequence); 
                    if ($self->give_qualloc) {
                        my $key;
                        my $qualloc = $detect->qualloc ($read_trim_static);
                        if($qualloc eq 'NA') {
                            $debug and warn 'Have an NA read' . $anchor_mate->name . "\n";
                        }
                        $key = $STRdetection2code_forw_strand{$qualloc}; # we have matched on the forward strand
                        if ($debug && $interesting_read) {
                            warn 'Found ok anchor on ' . ($anchor_direction == 1 ? 'forward' : 'reverse') . " strand, qualloc was $qualloc, with key $key\n";
                        }
                        # add the mate as a read in this category
                        if (defined ($key)) {
                            if (!exists($self->{$key})) { warn "No paired value for key $key." };
                            push @{$self->{$key}}, $detect_mate; 
                        }
                        if($debug && $interesting_read && !defined($key)) {
                            warn "The anchor read was found not to be suitable\n";
                        }
                    } 
                    if ($self->give_rep_in_read) {
                        my $rir = $detect->rep_in_read ($read_trim_static);
                        # give the inferred strand of the detect read, not the anchor
                        if(defined($rir)) {
                            my $rir_ob = STR::Rep_in_read->new (
                                a => $rir->[0],
                                b => $rir->[1],
                                c => $rir->[2],
                                strand => ($anchor_direction == 1 ? '-' : '+'),
                                bam_read => $detect_mate,
                            );
                            if (defined($rir)) {
                                push @{$self->rep_in_read}, $rir_ob;
                            }
                        }
                    } 
                }
            }
        }
    }
}

sub _find_mate {
    # find the mate of a given read
    my $read1 = $_[0];
    warn "  Finding mate for read " . $read1->name . "\n"; #TODO: verbose only
    my @reads2;
    if($read1->get_tag_values('M_UNMAPPED')) {
        # mate is unmapped, go find it!
        # something to do with $read1->name
        # Possible efficiency improvement: put all unmapped reads from original 
        #   region within the allowed boundaries into a hash and read from that.
        @reads2 = $read1->{sam}->features(
            -type => 'match',
            -name => $read1->name, # read should have the same name
            -seq_id => $read1->seq_id,
            -start => $read1->start,
            -end => $read1->start,
        );
    } else {
        # read is mapped somewhere far. Lets find it
        #  $read1->mate_start
        # and $read1->mate_seq_id
        @reads2 = $read1->{sam}->features(
            -name => $read1->name, # read should have the same name
            -seq_id => $read1->mate_seq_id,
            -start => $read1->mate_start,
            -end => $read1->mate_start,
        );
        # Don't accidently keep the original read (assuming it doesn't have the exact same start position)
        # If this same position occurs, we may have to look at the tag value directly
        @reads2 = grep { $_->start == $read1->mate_start } @reads2;
    }
    unless(@reads2 == 1) {
        if(@reads2 == 0) {
            warn "Could not find any mate for read " . $read1->name . "\n";
            return undef;
        } else {
            die "Multiple mates found for read " . $read1->name . "\n";
        }
    }
    return $reads2[0];
}

}

# class for the interpretation of CIGAR strings (probably not used)
## {
## package STR::Cigar;
## # Sometimes I just want to know the position that the end of the read is, and 
## # assume that soft clipped reads are entirely inside the last matching base.
## use Moo;
## extends 'Bio::Cigar';
## 
## sub qpos_to_rpos_soft {
##     my $self = shift;
##     my ($rpos, $op) = $self->qpos_to_rpos(@_);
##     if($op eq 'S') {
##         $rpos = $self->reference_length;
##     }
##     return wantarray ? ($rpos, $op) : $rpos;
## }
## 
## sub op_at_qpos_soft {
##     my $self = shift;
##     my (undef, $type) = $self->qpos_to_rpos_soft(@_);
##     return $type;
## }
## }


### 
{
package STR::Rep_in_read;
use Moose; 
use namespace::autoclean;
use autodie;

has [qw(a c)] => (
    isa => 'Natural0',
    is => 'rw',
);

has [qw(b)] => (
    isa => 'Natural',
    is => 'rw',
);

has 'strand' => (
    isa => 'Strand',
    is => 'rw',
);

has 'bam_read' => ( # keeps track of the read whose pair these values were generated from
    isa => 'Bio::DB::Bam::AlignWrapper',
    is => 'rw',
);

}

### 
{
package STR::Score;
use Moose; 
use namespace::autoclean;
use autodie;

has [qw(repeat_starts)] => (
    isa => 'Natural0',
    is => 'rw',
);

has [qw(mlength)] => ( # total number of bases we may start at
    isa => 'Natural',
    is => 'rw',
);

has 'bam_read' => ( # keeps track of the read whose pair these values were generated from
    isa => 'Bio::DB::Bam::AlignWrapper',
    is => 'rw',
);

}
1; # Exit success

