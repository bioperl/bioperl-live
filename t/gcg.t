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
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}

	use Test::More;
	plan tests => 5;
	use_ok('Bio::SeqIO');
	use_ok('Bio::Root::IO');
}

my $verbose = $ENV{'BIOPERLDEBUG'};

# test DOS linefeeds in gcg parser
my $str = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
								  ("t","data","test_badlf.gcg"),
								  -verbose => $verbose,
								  -format => 'GCG');
ok($str);
my $seq = $str->next_seq();
isa_ok ($seq, 'Bio::SeqI');
is(length($seq->seq), $seq->length);
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n" if( $verbose);
