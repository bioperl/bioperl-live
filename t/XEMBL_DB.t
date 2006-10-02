# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use constant NUMTESTS => 9;
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

	plan tests => NUMTESTS;

	unless( eval "require SOAP::Lite; require XML::DOM; 1;" ) {
      print STDERR "SOAP::Lite and/or XML::DOM not installed. This means that Bio::DB::XEMBL module is not usable. Skipping tests.\n";
      for( 1..NUMTESTS ) {
			skip("SOAP::Lite and/or XML::DOM not installed. This means that Bio::DB::XEMBL module is not usable. Skipping tests.\n",1);
      }
      $error = 1;
	}
}

if( $error ==  1 ) {
    exit(0);
}

END {
	foreach ( $Test::ntest..NUMTESTS) {
		skip('Cannot run XEMBL_DB tests',1);
	}
}

require Bio::DB::XEMBL;

my $testnum;
my $verbose = 1;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my ($db,$seq,$seqio);
# get a single seq

$seq = $seqio = undef;
eval {
ok defined($db = new Bio::DB::XEMBL(-verbose=>$verbose)); 
ok(defined($seq = $db->get_Seq_by_acc('J00522')));
ok( $seq->length, 408);
ok(defined($seq = $db->get_Seq_by_acc('J02231')));
ok $seq->id, 'BUM';
ok( $seq->length, 200); 
ok(defined($seqio = $db->get_Stream_by_batch(['BUM'])));
undef $db; # testing to see if we can remove gb
ok( defined($seq = $seqio->next_seq()));
ok( $seq->length, 200);
};
if( $@ ) { 
  skip('Skip server may be down',1);
  exit(0);
}

$seq = $seqio = undef;

eval {
    $db = new Bio::DB::XEMBL(-verbose => $verbose,
			    -retrievaltype => 'tempfile',
			    -format => 'fasta'
			    ); 
    ok( defined($seqio = $db->get_Stream_by_batch(['J00522 AF303112 J02231'])));
    ok($seqio->next_seq->length, 408);
    ok($seqio->next_seq->length, 1611);
    ok($seqio->next_seq->length, 200);
};

if ($@) {
    warn "Batch access test failed.\nError: $@\n";
    foreach ( $Test::ntest..NUMTESTS ) { skip('no network access',1); }
}
