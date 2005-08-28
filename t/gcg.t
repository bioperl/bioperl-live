# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}

	use Test;
	plan tests => 4;
}

if( $error == 1 ) {
	exit(0);
}

use Bio::SeqIO;
use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

# test DOS linefeeds in gcg parser
my $str = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
								  ("t","data","test_badlf.gcg"),
								  -verbose => $verbose,
								  -format => 'GCG');
ok($str);
ok ( my $seq = $str->next_seq());
ok(length($seq->seq) > 0 );
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n" if( $verbose);
