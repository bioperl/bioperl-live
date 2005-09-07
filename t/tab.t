# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
BEGIN {
	$NUMTESTS = 7;
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

use Bio::SeqIO;
use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $io = Bio::SeqIO->new(-format => 'tab',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data test.tab) ));

while (my $seq = $io->next_seq) {
	ok defined $seq;
	ok($seq->length, 358);
	ok($seq->display_id =~ /^roa\d_drome$/);
}
