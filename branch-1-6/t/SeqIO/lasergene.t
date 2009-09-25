# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 13);
	
	use_ok('Bio::SeqIO::lasergene');
}

my $verbose = test_debug();

#
# Positive tests
#

my $io = Bio::SeqIO->new(
  -format => 'lasergene',
  -verbose => $verbose,
  -file => test_input_file('test.lasergene')
);

ok($io);
isa_ok($io, 'Bio::SeqIO');

my $seq;

ok($seq = $io->next_seq);
is($seq->length, 12*3);
is($seq->subseq(1,12), 'ATCGATCGATCG');

ok($seq = $io->next_seq);
is($seq->length, 200);

ok($seq = $io->next_seq);
is($seq->length, 70*5+12);

ok(not defined $io->next_seq);

#
# Negative tests
#

$io = Bio::SeqIO->new(
  -format => 'lasergene',
  -verbose => $verbose,
  -file => test_input_file('test.fasta') # not lasergene!
);

ok($io);

eval { 
  $io->next_seq;
};
like($@, qr/unexpected end of file/i);
