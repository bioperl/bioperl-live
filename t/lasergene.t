# -*-Perl-*-
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
BEGIN {
	$NUMTESTS = 13;
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test::More; };
	if ( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => $NUMTESTS;
	use_ok('Bio::SeqIO::raw');
	use_ok('Bio::Root::IO');	
}

my $verbose = $ENV{'BIOPERLDEBUG'};

#
# Positive tests
#

my $io = Bio::SeqIO->new(
  -format => 'lasergene',
  -verbose => $verbose,
  -file => Bio::Root::IO->catfile(qw(t data test.lasergene))
);

ok($io);

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
  -file => Bio::Root::IO->catfile(qw(t data test.fasta)) # not lasergene!
);

ok($io);

eval { 
  $io->next_seq;
};
like($@, qr/unexpected end of file/i);

