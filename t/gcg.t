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

# test DOS linefeeds in gcg parser
my $str = Bio::SeqIO->new(-file => test_input_file('test_badlf.gcg'),
								  -verbose => $verbose,
								  -format => 'GCG');
ok($str);
my $seq = $str->next_seq();
isa_ok ($seq, 'Bio::SeqI');
is(length($seq->seq), $seq->length);
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n" if( $verbose);
