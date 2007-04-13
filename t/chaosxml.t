# -*-Perl-*-
# $Id$
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($error $NUMTESTS $out_file);
BEGIN {
	$NUMTESTS = 3;
	use File::Spec;
	$out_file = File::Spec->catfile(qw(t data tmp-chaosxml));
	eval { require Test::More; };
	if ( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	eval {
		require Data::Stag;
	};
	if ( $@ ) {
		plan skip_all => "Data::Stag::XMLWriter not installed, cannot perform chaosxml tests";
	} else {
		plan tests => $NUMTESTS;
	}
	use_ok('Bio::SeqIO');
	use_ok('Bio::Root::IO');
}

END {
	unlink $out_file if -e $out_file;
}

my $verbose = $ENV{'BIOPERLDEBUG'};

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
