# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
   use lib '.';
   use Bio::Root::Test;
   
   test_begin(-tests => 260);
   
   use_ok('Bio::SeqIO::genbank');
}

my $verbose = test_debug();

my $ast = Bio::SeqIO->new(-format  => 'gbxml',
                          -verbose => $verbose,
                          -file    => test_input_file('roa1.gbxml'));
isa_ok($ast, 'Bio::SeqIO');
$ast->verbose($verbose);
my $as = $ast->next_seq();
is $as->molecule, 'mRNA',$as->accession_number;
is $as->alphabet, 'dna';
is($as->primary_id, 3598416);
my @class = $as->species->classification;
is $class[$#class],'Eukaryota';


