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

    $NUMTESTS = 23;
    plan tests => $NUMTESTS;

    eval {
	require IO::String; 
	require LWP::UserAgent;
	require Bio::WebAgent;
    }; 
    if( $@ ) {
        warn("IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests");
	$ERROR = 1;
    }
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
	skip('unable to run all of the tests depending on web access',1);
    }
}

exit 0 if $ERROR ==  1;

use Data::Dumper;

require Bio::Tools::Analysis::DNA::ESEfinder;
use Bio::PrimarySeq;
use Bio::Seq;
require Bio::WebAgent;

ok 1;

my $verbose = 0;
$verbose = 1 if $DEBUG;
my $tool;


#######all these tests work with 1ary seq########
my $seq = Bio::PrimarySeq->new(-id=>'bioperl',
                               -seq=>'atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt');

ok $tool = Bio::Tools::Analysis::DNA::ESEfinder->new(-verbose =>$verbose, -seq => $seq);
if ( $DEBUG ) {
    ok $tool->run ( );
    
    ok my @res = $tool->result('Bio::SeqFeatureI');
    #new tests her in v 1.2
    ok my $raw = $tool->result('');
    ok my $parsed = $tool->result('parsed');
    ok my $meta = $tool->result('all');
    ok ($parsed->[0][1], 41);
    if (scalar @res > 0) {
	ok 1;
    } else {
	skip('No network access - could not connect to ESEfinder server', 1);
    }
    if (!$METAERROR) { #if Bio::Seq::Meta::Array available
	ok($meta->{'seq'}, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt" );
	ok( $meta->named_submeta_text('ESEfinder_SRp55', 1,2), "-3.221149 -1.602223");
	ok ($meta->seq, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt" );
    }
######## now repeat with Bio::Seq object, metasequence tests fail ########
    $seq = Bio::Seq->new(-id=>'bioperl',
			 -seq=>'atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt');
    
    ok $tool = Bio::Tools::Analysis::DNA::ESEfinder->new(-verbose =>$verbose, -seq => $seq);
    ok $tool->run ( );
    
    ok @res = $tool->result('Bio::SeqFeatureI');
#new tests her in v 1.2
    ok $raw = $tool->result('');
    ok $parsed = $tool->result('parsed');
    ok $meta = $tool->result('all');
    ok ($parsed->[0][1], 41);
    if (scalar @res > 0) {
	ok 1;
    } else {
	skip('No network access - could not connect to ESEfinder server', 1);
    }

    if (!$METAERROR) { #if Bio::Seq::Meta::Array available
	skip("meta sequence returns undef with Bio::Seq object ",$meta->{'seq'}, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt" );
	skip( "meta sequence returns undef with Bio::Seq object ", $meta->named_submeta_text('ESEfinder_SRp55', 1,2), "-3.221149 -1.602223");
	skip ("meta sequence returns undef with Bio::Seq object ", $meta->seq, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt" );
    }
} else { 
    foreach ( $Test::ntest..$NUMTESTS) {
	skip('skipping tests to avoid timeouts - set BIOPERLDEBUG env variable to 1 to run.',1);
    }
}
