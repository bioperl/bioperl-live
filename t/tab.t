# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 7);
	
	use_ok('Bio::SeqIO');
}

my $verbose = test_debug();

my $io = Bio::SeqIO->new(-format => 'tab',
								 -verbose => $verbose,
								 -file => test_input_file('test.tab'));

while (my $seq = $io->next_seq) {
	ok ( $seq && defined $seq, 'seq is defined' ) ;
	is ( $seq->length, 358, 'check seq length'  ) ;
	like ($seq->display_id, qr/^roa\d_drome$/, 'check matching' );
}
