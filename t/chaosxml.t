# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($error $NUMTESTS $out_file);
BEGIN {
	$NUMTESTS = 2;
	$error = 0;
	use File::Spec;
	$out_file = File::Spec->catfile(qw(t data tmp-chaosxml));
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
		warn "Data::Stag::Writer not installed, cannot perform chaosxml tests\n";
   } 
	use Test;
	plan tests => $NUMTESTS;
}

END { 
   foreach ( $Test::ntest..$NUMTESTS) {
      skip('Unable to run all of the chaosxml tests',1);
   }
}


if ( $error == 1 ) {
  exit(0);
}

use Bio::SeqIO;
use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

# currently chaosxml is write-only
my $in = Bio::SeqIO->new(-format => 'genbank',
								 -verbose => $verbose,
								 -file => Bio::Root::IO->catfile
								 qw(t data AE003644_Adh-genomic.gb) );

my $seq = $in->next_seq;

my $out = Bio::SeqIO->new(-file => ">$out_file",
								  -verbose => $verbose,
								  -format => 'chaosxml');
$out->write_seq($seq);
ok (-e $out_file);

END {
	unlink $out_file if -e $out_file;
}
