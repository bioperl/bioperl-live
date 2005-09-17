# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: Sopma.t,v 1.1 2003/07/23 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use strict;
use vars qw($NUMTESTS $DEBUG $ERROR $METAERROR);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	$ERROR = 0;
	if( $@ ) {
		use lib 't';
	}
	use Test;

	$NUMTESTS = 15;
	plan tests => $NUMTESTS;

	eval {
		require IO::String; 
		require LWP::UserAgent;
	};
	if( $@ ) {
		warn("IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests");
	$ERROR = 1;
	}
	# check this is available, set error flag if not.
	eval {
		require Bio::Seq::Meta::Array;
	};
	if ($@) {
		warn ("Bio::Seq::Meta::Array not installed - will skip tests using meta sequences");
		$METAERROR = 1;
	}
}

END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('unable to run all of the Sopma tests',1);
	}
}

exit 0 if $ERROR ==  1;

use Data::Dumper;
use Bio::PrimarySeq;
require Bio::Tools::Analysis::Protein::Sopma;

ok 1;

my $verbose = 0;
$verbose = 1 if $DEBUG;

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

my $seq = Bio::PrimarySeq->new(
  -seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
  -display_id => 'test2');
ok $tool = Bio::Tools::Analysis::Protein::Sopma->new( -seq=>$seq,
                                                      -window_width => 15);
if( $DEBUG ) {
	ok $tool->run ();
	exit if $tool->status eq 'TERMINATED_BY_ERROR';
	ok my $raw = $tool->result('');
	ok my $parsed = $tool->result('parsed');
	ok ($parsed->[0]{'helix'}, '104');
	ok my @res = $tool->result('Bio::SeqFeatureI');
	ok my $meta = $tool->result('meta', "ww15");

	ok $tool->window_width(21);
	ok $tool->clear();
	ok $tool->run;
	ok my $meta2 = $tool->result('meta', "ww21");
	if (!$METAERROR) { 
		# if Bio::Seq::Meta::Array available
		# meta sequence contains data...
		# but not available thru method call...??
		ok ($meta->named_submeta_text('Sopma_helix|ww15',1,2), '104 195');
		ok ($meta->seq, 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS');
	}
} else {
	for ( $Test::ntest..$NUMTESTS) {
		skip("Skipping tests which require remote servers - set env variable BIOPERLDEBUG to test",1);
	}
}
