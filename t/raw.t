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
	eval { require Test::More ; };
	if ( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => $NUMTESTS;
}

use_ok('Bio::SeqIO::raw','Bio::SeqIO::raw can be used' ) ;
use_ok('Bio::Root::IO', 'Bio::Root::IO can be used' ) ;

my $verbose = $ENV{'BIOPERLDEBUG'};

my $io = Bio::SeqIO->new(-format => 'raw',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data test.raw) ));

my $seq ; 

ok( $seq = $io->next_seq, 'got next seq');
is( $seq->length, 358, 'checked seq length');
ok( $seq = $io->next_seq, 'got next seq');
is( $seq->length, 158, 'checked seq length');
