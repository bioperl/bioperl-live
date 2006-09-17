# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$ 

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

    $NUMTESTS = 12;
    plan tests => $NUMTESTS;

    eval {
	require IO::String; 
	require LWP::UserAgent;
	require Bio::WebAgent;
	require HTML::HeadParser;
	require HTTP::Request::Common;
	require Bio::Tools::Analysis::DNA::ESEfinder;
        1;
    }; 
    if( $@ ) {
        warn("IO::String, LWP::UserAgent, Bio::WebAgent, HTML::HeadParser, or HTTP::Request::Common not installed. This means that the module is not usable. Skipping tests\n");
	$ERROR = 1;
    } else {
     eval {
	require Bio::Seq::Meta::Array; 1;
     };
     if ($@) {
	warn ("Bio::Seq::Meta::Array not installed - will skip tests using meta sequences");
	$METAERROR = 1;
    }
  }
}

END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('unable to complete ESEfinder tests',1);
	}
}

if($ERROR ==  1 ) {
  warn("exiting early\n");
  exit(0);
}

use Data::Dumper;
use Bio::PrimarySeq;
use Bio::Seq;

ok 1;

my $verbose = $DEBUG;
my $tool;


#######all these tests work with 1ary seq########
my $seq = Bio::PrimarySeq->new(-id=>'bioperl',
                               -seq=>'atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt');

ok $tool = Bio::Tools::Analysis::DNA::ESEfinder->new(-verbose =>$verbose, -seq => $seq);
if( $DEBUG ) {
	eval {
		ok $tool->run;
	};
	if ($@) {
		foreach ( $Test::ntest..$NUMTESTS) { 
			skip('Could not connect to ESEfinder server', 1);
		}
		exit(0);
	}
	
    ok my @res = $tool->result('Bio::SeqFeatureI');
	ok @res > 0;
    ok my $raw = $tool->result('');
    ok my $parsed = $tool->result('parsed');
    ok my $meta = $tool->result('all');
    ok ($parsed->[0][1], 41);
	
    if (!$METAERROR) {              #if Bio::Seq::Meta::Array available
	ok($meta->{'seq'}, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt" );
	ok( $meta->named_submeta_text('ESEfinder_SRp55', 1,2), "-3.221149 -1.602223");
	ok ($meta->seq, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt" );
    }
} else { 
   for ( $Test::ntest..$NUMTESTS) {
	skip("Skipping tests which require remote servers - set env variable BIOPERLDEBUG to test",1);
    }
}
