# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
BEGIN {
	$NUMTESTS = 5;
	eval { require Test::More; };
	if ( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => $NUMTESTS;
	use_ok('Bio::SeqIO::metafasta');
	use_ok('Bio::Root::IO');
}

my $verbose = $ENV{'BIOPERLDEBUG'};

my $io = Bio::SeqIO->new(-format => 'metafasta',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data test.metafasta) ));

ok(my $seq = $io->next_seq);
is($seq->seq, "ABCDEFHIJKLMNOPQRSTUVWXYZ");
is($seq->display_id,'test');
