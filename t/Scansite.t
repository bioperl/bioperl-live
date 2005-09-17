# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: Scansite.t,v 1.1 2003/11/18 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG $ERROR);
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

	$NUMTESTS = 12;
	plan tests => $NUMTESTS;

	eval {
		require IO::String; 
		require LWP::UserAgent;
	};
	if( $@ ) {
		warn("IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests");
		$ERROR = 1;
	}
}

END {
	if ($DEBUG) {
		foreach ( $Test::ntest..$NUMTESTS) {
			skip('unable to run all of the Scansite tests',1);
		}
	} else {
		foreach ( $Test::ntest..$NUMTESTS) {
			skip('set env BIOPERLDEBUG to run tests over web',1);
		}
	}
}

exit 0 if $ERROR ==  1;

use Data::Dumper;

require Bio::Tools::Analysis::Protein::Scansite;
use Bio::SeqIO;
use Bio::PrimarySeq;
require Bio::WebAgent;

ok 1;

my $verbose = 0;
$verbose = 1 if $DEBUG;

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);


my $seqio=new Bio::SeqIO( -verbose => $verbose,
                  -format => 'swiss',
                  -file   => Bio::Root::IO->catfile('t','data', 'swiss.dat'));

my $seq = $seqio->next_seq();
ok $tool = Bio::Tools::Analysis::Protein::Scansite->new( 
					-seq=>$seq->primary_seq);
ok $tool->stringency('Low');
ok $tool->stringency(), 'Low';
ok $tool->protein_id(), $tool->seq->display_id();
exit unless $DEBUG;
ok $tool->run ();
exit if $tool->status eq 'TERMINATED_BY_ERROR';
ok my $raw = $tool->result('');
print $raw if $verbose;
ok my $parsed = $tool->result('parsed');
ok $parsed->[0]{'site'}, 'T101';
ok my @res = $tool->result('Bio::SeqFeatureI');
ok $res[0]->start, 101;
