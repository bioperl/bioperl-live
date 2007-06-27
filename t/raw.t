# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 5);
	
	use_ok('Bio::SeqIO');
}

my $verbose = test_debug();

my $io = Bio::SeqIO->new(-format => 'raw',
								 -verbose => $verbose,
								 -file => test_input_file('test.raw'));

my $seq ; 

ok( $seq = $io->next_seq, 'got next seq');
is( $seq->length, 358, 'checked seq length');
ok( $seq = $io->next_seq, 'got next seq');
is( $seq->length, 158, 'checked seq length');
