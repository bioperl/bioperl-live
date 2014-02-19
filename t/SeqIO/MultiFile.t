# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 5);

    use_ok 'Bio::SeqIO::MultiFile';
}

my $verbose = test_debug();


# Test multiple files, with a specified format
ok my $mf = Bio::SeqIO::MultiFile->new(
    -format  => 'Fasta' ,
    -verbose => $verbose,
    -files   => [ test_input_file('multi_1.fa'), test_input_file('multi_2.fa')],
);

my $count = 0;
while (my $seq = $mf->next_seq() ) {
    $count++;
}
is $count, 12;


# Automatically determine format
ok $mf = Bio::SeqIO::MultiFile->new(
    -verbose => $verbose,
    -files   => [ test_input_file('multi_1.fa'), test_input_file('multi_2.fa')],
);

is $mf->format, 'fasta';
