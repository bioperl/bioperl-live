# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 6);
	
	use_ok('Bio::SeqIO');
}

# There are many other tests of fasta I/O in t/, this
# is a dedicated script that could be further customized
 
my $verbose = test_debug();

my $io = Bio::SeqIO->new(-format => '',
								 -verbose => $verbose,
								 -file => test_input_file('test.fasta'));

ok(my $seq = $io->next_seq);
is($seq->length, 358);
is($seq->display_id,'roa1_drome');
is($seq->desc,'Rea guano receptor type III >> 0.1');
is($seq->alphabet,'protein');
