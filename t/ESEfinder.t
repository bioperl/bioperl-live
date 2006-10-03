# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$ 

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);

BEGIN {
	$NUMTESTS = 15;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	
	eval {require Test::More;};
	if ($@) {
		use lib 't';
	}
	use Test::More;
	
	eval {
		require IO::String; 
		require LWP::UserAgent;
		require Bio::WebAgent;
		require HTML::HeadParser;
		require HTTP::Request::Common;
	};
	if ($@) {
		plan skip_all => 'IO::String, LWP::UserAgent, Bio::WebAgent, HTML::HeadParser, or HTTP::Request::Common not installed. This means that the module is not usable. Skipping tests';
	}
	else {
		plan tests => $NUMTESTS;
	}
	
	use_ok('Bio::Tools::Analysis::DNA::ESEfinder');
	use_ok('Data::Dumper');
	use_ok('Bio::PrimarySeq');
    use_ok('Bio::Seq');
}

#######all these tests work with 1ary seq########
ok my $seq = Bio::PrimarySeq->new(-id=>'bioperl',
                               -seq=>'atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt');
ok my $tool = Bio::Tools::Analysis::DNA::ESEfinder->new(-seq => $seq);

SKIP: {
	skip "Skipping tests which require remote servers, set BIOPERLDEBUG=1 to test", 9 unless $DEBUG;
	eval {ok $tool->run;};
	skip "Could not connect to ESEfinder server, skipping those tests", 9 if $@;
    ok my @res = $tool->result('Bio::SeqFeatureI');
    ok @res > 0;
    ok my $raw = $tool->result('');
    ok my $parsed = $tool->result('parsed');
    ok my $meta = $tool->result('all');
    is $parsed->[0][1], 41;
	
    eval {require Bio::Seq::Meta::Array;};
	skip "Bio::Seq::Meta::Array not installed. Skipping tests using meta sequences", 3 if $@;
	is $meta->{'seq'}, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt";
	is $meta->named_submeta_text('ESEfinder_SRp55', 1,2), "-3.221149 -1.602223";
	is $meta->seq, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt";
}
