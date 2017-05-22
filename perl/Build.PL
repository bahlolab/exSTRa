use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Bio::STR::exSTRa',
    license             => 'perl',
    dist_abstract       => 'Prepare expanded repeat data from BAM files for the exSTRa R package',
    dist_author         => 'Rick Tankard <tankard@wehi.edu.au>',
    build_requires => {
        'Test::More' => '0.10',
    },
);

$builder->create_build_script();

