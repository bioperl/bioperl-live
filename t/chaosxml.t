# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($error $NUMTESTS);
BEGIN {
	$NUMTESTS = 2;
	$error = 0;
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	#
	eval {
		require Data::Stag;
	};
	if ( $@ ) {
		$error = 1;
		warn "Data::Stag not installed, cannot perform chaosxml tests\n";
   } 
	use Test;
	plan tests => $NUMTESTS;
}

if ( $error == 1 ) {
  exit(0);
}

use Bio::SeqIO;
use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $io = Bio::SeqIO->new(-format => 'chaosxml',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data Rab1.chaos-xml) ));
ok(my $seq = $io->next_seq);
# ok($seq->length, 1063);
# ok($seq->display_id,'AE00373');
# ok($seq->molecule,'dna');

END { 
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Unable to run all of the chaosxml tests',1);
   }
}
