# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
use lib '.','./blib/lib';

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

    $NUMTESTS = 9;
    plan tests => $NUMTESTS;

    unless( eval "use IO::String; 1;" ) {
      warn $@;
      print STDERR "IO::String not installed. This means that Bio::DB::BioFetch module is not usable. Skipping tests.\n";
      for( 1..$NUMTESTS ) {
	skip(1,"IO::String not installed. This means that Bio::DB::BioFetch module is not usable. Skipping tests.\n");
      }
      $error = 1;
    }
}

if( $error ==  1 ) {
    exit(0);
}

require Bio::DB::BioFetch;

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my ($db,$seq,$seqio);
# get a single seq

$seq = $seqio = undef;

ok defined($db = new Bio::DB::BioFetch);
ok(defined($seq = $db->get_Seq_by_acc('J00522')));
ok( $seq->length, 408);
ok(defined($seq = $db->get_Seq_by_acc('J02231')));
ok $seq->id, 'BUM';
ok( $seq->length, 200); 
ok(defined($seqio = $db->get_Stream_by_id(['BUM'])));
undef $db; # testing to see if we can remove gb
ok( defined($seq = $seqio->next_seq()));
ok( $seq->length, 200);

exit;
$seq = $seqio = undef;

eval {
    $db = new Bio::DB::BioFetch(-retrievaltype => 'tempfile',
				-format => 'fasta'
			       );
    ok( defined($seqio = $db->get_Stream_by_batch(['J00522 AF303112 J02231'])));
    ok($seqio->next_seq->length, 408);
    ok($seqio->next_seq->length, 1611);
    ok($seqio->next_seq->length, 200);
};

if ($@) {
    warn "Batch access test failed.\nError: $@\n";
    foreach ( $Test::ntest..$NUMTESTS ) { skip('no network access',1); }
}


