# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
BEGIN {
	$NUMTESTS = 5;
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	use Test;
	plan tests => $NUMTESTS;
}

use Bio::SeqIO::raw;
use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $io = Bio::SeqIO->new(-format => 'raw',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data test.raw) ));

ok(my $seq = $io->next_seq);
ok($seq->length, 358);
ok($seq = $io->next_seq);
ok($seq->length, 158);
