# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
my $error;

BEGIN { 
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	$error = 0;
	if( $@ ) {
		use lib 't';
	}
	use Test;

	$NUMTESTS = 13;
	plan tests => $NUMTESTS;
	eval { require IO::String; 
			 require LWP::UserAgent; 
			 require HTTP::Request::Common;
		 };
	if( $@ ) {
		for( $Test::ntest..$NUMTESTS ) {
			skip("IO::String,LWP::UserAgent, or HTTP::Request::Common not installed. This means the Bio::DB::* modules are not usable. Skipping tests",1);
		}
		$error = 1;
	}
}

END {
	for ( $Test::ntest..$NUMTESTS ) {
		skip("Unable to complete RefSeq tests - set env variable BIOPERLDEBUG to test",1);
	}
}

if( $error ==  1 ) {
    exit(0);
}

require  Bio::DB::RefSeq;
require  Bio::DB::GenBank;
require  Bio::DB::EMBL;

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my ($db,$seq,$db2,$seq2,$seqio);
# get a single seq

$seq = $seqio = undef;

#test redirection from GenBank and EMBL
$verbose = -1;
#GenBank
ok $db = Bio::DB::GenBank->new('-verbose'=>$verbose) ;     
#EMBL
ok $db2 = Bio::DB::EMBL->new('-verbose'=>$verbose) ;     

eval {
    $seq = $db->get_Seq_by_acc('NT_006732');
    $seq2 = $db2->get_Seq_by_acc('NT_006732');
};
ok $@;

exit unless $DEBUG;
eval {
    ok($seq = $db->get_Seq_by_acc('NM_006732'));
    ok($seq && $seq->length eq 3775);
    ok $seq2 = $db2->get_Seq_by_acc('NM_006732');
    ok($seq2 && $seq2->length eq 3775);
};

if ($@) {
    if( $DEBUG ) {
	print STDERR "Warning: Couldn't connect to RefSeq with Bio::DB::RefSeq.pm!\n" . $@;
    }
    exit(0);
}



$verbose = 0;

eval { 
    ok defined($db = Bio::DB::RefSeq->new(-verbose=>$verbose)); 
    ok(defined($seq = $db->get_Seq_by_acc('NM_006732')));
    ok( $seq->length, 3775);
    ok defined ($db->request_format('fasta'));
    ok(defined($seq = $db->get_Seq_by_acc('NM_006732')));
    ok( $seq->length, 3775); 
};

if ($@) {
    if( $DEBUG ) {
	print STDERR "Warning: Couldn't connect to RefSeq with Bio::DB::RefSeq.pm!\n" . $@;
    }
    exit(0);
}

