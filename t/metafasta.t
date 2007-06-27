# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 4);
	
	use_ok('Bio::SeqIO');
}

my $verbose = test_debug();

my $io = Bio::SeqIO->new(-format => 'metafasta',
								 -verbose => $verbose,
								 -file => test_input_file('test.metafasta'));

ok(my $seq = $io->next_seq);
is($seq->seq, "ABCDEFHIJKLMNOPQRSTUVWXYZ");
is($seq->display_id,'test');
