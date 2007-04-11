# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);

BEGIN {

	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test::More; };
	
	if ( $@ ) {
		use lib 't/lib';
	}
    
    $NUMTESTS = 8 ;
	use Test::More ;
	plan tests => $NUMTESTS ;
	
}

use_ok ('Bio::SeqIO', 'Bio::SeqIO can be used') ;
use_ok('Bio::Root::IO', 'Bio::Root::IO can be used') ;

my $verbose = $ENV{'BIOPERLDEBUG'};

my $io = Bio::SeqIO->new(-format => 'tab',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data test.tab) ));

while (my $seq = $io->next_seq) {
	ok ( $seq && defined $seq, 'seq is defined' ) ;
	is ( $seq->length, 358, 'check seq length'  ) ;
	like ($seq->display_id, qr/^roa\d_drome$/, 'check matching' );
}
