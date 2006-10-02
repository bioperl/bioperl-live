# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# This test file is adapted from EMBL_DB.t

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl SeqHound_DB.t'

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

my $error;

BEGIN {
	use Test;

	$NUMTESTS = 14;
	plan tests => $NUMTESTS;
	eval { require IO::String;
			 require LWP::UserAgent;
		 };
	if( $@ ) {
		for( $Test::ntest..$NUMTESTS ) {
			skip("IO::String or LWP::UserAgent not installed. This means the Bio::DB::* modules are not usable. Skipping tests",1);
		}
		$error = 1;
	}
}
END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Unable to run all of the SeqHound_DB tests',1);
	}
	unlink $Bio::DB::SeqHound::LOGFILENAME if -f $Bio::DB::SeqHound::LOGFILENAME;
}

if( $error ) {
    exit(0);
}

use Bio::DB::SeqHound;

my $verbose = -1;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my ($db,$seq,$seqio);
# get a single seq

$seq = $seqio = undef;

eval { 
    $db = new Bio::DB::SeqHound(-verbose=>$verbose);
    ok(defined($seq = $db->get_Seq_by_acc('J00522')));
    ok( $seq->length, 408); 
    ok defined ($db->request_format('fasta'));
    ok(defined($seq = $db->get_Seq_by_acc('NP_862707')));
    ok( $seq->accession, 'NP_862707');
    ok( $seq->length, 227); 
    ok( defined($db = new Bio::DB::SeqHound(-verbose=>$verbose, 
					-retrievaltype => 'tempfile')));
    ok(defined($seqio = $db->get_Stream_by_id(['BTACHRE'])));
    undef $db; # testing to see if we can remove db
    ok( defined($seq = $seqio->next_seq()));
    ok( $seq->length, 1621);
};

if ($@) {
    if( $DEBUG ) {
	warn "Warning: Couldn't connect to SeqHound with Bio::DB::SeqHound.pm!\n$@";
    }
    foreach ( $Test::ntest..$NUMTESTS) { 
	 skip('Skip could not connect to seqhound',1);
     }
    exit(0);
}

$seq = $seqio = undef;

eval {
    $db = Bio::DB::SeqHound->new(-verbose => $verbose,
			    -retrievaltype => 'tempfile',
			    -format => 'genbank'
			    ); 
    ok( defined($seqio = $db->get_Stream_by_acc(['J00522', 'AF303112', 'J02231'])));
	ok($seqio->next_seq->length, 408);
    ok($seqio->next_seq->length, 1611);
    ok($seqio->next_seq->length, 200);
};

if ($@) {
    if( $DEBUG ) {
	warn "Batch access test failed.\nError: $@\n";
    }
    foreach ( $Test::ntest..$NUMTESTS ) { skip('could not connect to seqhound',1); }
}


