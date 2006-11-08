# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$ 
use strict;
use vars qw($NUMTESTS $DEBUG);

BEGIN {
	$NUMTESTS = 14;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	
	eval {require Test::More;};
	if ($@) {
		use lib 't/lib';
	}
	use Test::More;
	
	eval {
		require IO::String; 
		require LWP::UserAgent;
	};
	if ($@) {
		plan skip_all => 'IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests';
	}
	else {
		plan tests => $NUMTESTS;
	}
	
	use_ok('Bio::Tools::Analysis::Protein::NetPhos');
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::WebAgent');
}

my $verbose = 0;
$verbose = 1 if $DEBUG;

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

SKIP: {
	skip "Skipping tests which require network access, set BIOPERLDEBUG=1 to test", 10 unless $DEBUG;
	ok $tool->sleep;
	is $tool->delay(1), 1;
	ok $tool->sleep;
	ok $tool->timeout(120); # LWP::UserAgent method
	is $tool->url('http://a.b.c/'), 'http://a.b.c/';
	
	
	my $seq = Bio::PrimarySeq->new(-id=>'bioperl',
								   -seq=>'ABCDEFGHIJKLLKJFHSAKNDJFPSINCSJNDSKNSN');
	
	ok $tool = Bio::Tools::Analysis::Protein::NetPhos->new(-verbose =>$verbose);
	$tool->timeout(15);
	ok $tool->run ( {seq=>$seq, threshold=>0.9} );
	if ($tool->status eq 'TERMINATED_BY_ERROR') {
		skip "Running of the tool was terminated by an error, probably network/ NetPhos server error", 3;
	}
	my @res = $tool->result('Bio::SeqFeatureI');
	unless (@res) {
		skip "Didn't get any results from NetPhos server, probable network/server error", 3;
	}
	#new tests her in v 1.2
	ok my $raw = $tool->result('');
	ok my $parsed = $tool->result('parsed');
	is $parsed->[0][1], '0.934';
}
