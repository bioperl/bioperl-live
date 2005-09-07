# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
BEGIN {
	$NUMTESTS = 6;
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	#
	use Test;
	plan tests => $NUMTESTS;
}

use Bio::SeqIO::fasta;
use Bio::Root::IO;

# There are many other tests of fasta I/O in t/, this
# is a dedicated script that could be further customized
 
my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $io = Bio::SeqIO->new(-format => '',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data test.fasta) ));

ok(my $seq = $io->next_seq);
ok($seq->length, 358);
ok($seq->display_id,'roa1_drome');
ok($seq->desc,'Rea guano receptor type III >> 0.1');
ok($seq->alphabet,'protein');
