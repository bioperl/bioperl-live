# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
BEGIN {
	$NUMTESTS = 8;
	eval { require Test::More; };
	if ( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => $NUMTESTS;
	use_ok('Bio::SeqIO::fasta');
	use_ok('Bio::Root::IO');
}

# There are many other tests of fasta I/O in t/, this
# is a dedicated script that could be further customized
 
my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $io = Bio::SeqIO->new(-format => '',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 (qw(t data test.fasta) ));

ok(my $seq = $io->next_seq);
is($seq->length, 358);
is($seq->display_id,'roa1_drome');
is($seq->desc,'Rea guano receptor type III >> 0.1');
is($seq->alphabet,'protein');
