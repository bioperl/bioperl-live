# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);

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

    $NUMTESTS = 51;
    plan tests => $NUMTESTS;
    eval { require 'IO/String.pm' };
    if( $@ ) {
	print STDERR "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n";
	for( 1..$NUMTESTS ) {
	    skip(1,"IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests");
	}
       $error = 1; 
    }
}

if( $error ==  1 ) {
    exit(0);
}

require Bio::DB::GenBank;
require Bio::DB::GenPept;
require Bio::DB::SwissProt;

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


my ($gb,$seq,$seqio);
# get a single seq



eval {         
    ok defined ( $gb = new Bio::DB::GenBank('-verbose'=>$verbose) );     
    ok( defined ($seq = $gb->get_Seq_by_id('MUSIGHBA1')));
    ok( $seq->length, 408); 
    ok( defined ($seq = $gb->get_Seq_by_acc('AF303112'))); 
    ok($seq->length, 1611);
    ok( defined ($seq = $gb->get_Seq_by_version('AF303112.1'))); 
    ok($seq->length, 1611);
    ok( defined ($seq = $gb->get_Seq_by_gi('405830'))); 
    ok($seq->length, 1743);
};
if ($@) {
    warn "Warning: Couldn't connect to Genbank with Bio::DB::GenBank.pm!\nError: $@\nDo you have network access? Skipping all other tests";
    foreach ( $Test::ntest..$NUMTESTS ) { skip(1,1, 'no network access'); }
    exit(0);
}

$seq = $seqio = undef;

eval {
    ok( defined($seqio = $gb->get_Stream_by_batch([ qw(J00522 AF303112 
							 2981014)])));
    ok($seqio->next_seq->length, 408);
    ok($seqio->next_seq->length, 1611);
    ok($seqio->next_seq->length, 1156);
};

if ($@) {
    warn "Batch access test failed.\nError: $@\n";
    foreach ( $Test::ntest..$NUMTESTS ) { skip(1,1,'no network access'); }
    exit(0);
}
$seq = $seqio = undef;

eval { 
    ok defined($gb = new Bio::DB::GenPept('-verbose'=>$verbose)); 
    ok( defined($seq = $gb->get_Seq_by_id('195055')));
    ok( $seq->length, 136); 
    $seq = $gb->get_Seq_by_acc('AAC06201');
    ok(defined $seq);
    ok($seq->length, 353);
    $seqio = $gb->get_Stream_by_batch([ qw(AAC06201 195055)]);
    ok( defined $seqio);
    ok( $seqio->next_seq->length(), 353);
    ok( $seqio->next_seq->length(), 136);
};

if ($@) {
    warn "Warning: Couldn't connect to Genbank with Bio::DB::GenPept.pm!\n$@";
    foreach( $Test::ntest..$NUMTESTS ) { 
	skip(1,1,1,'could not connect with GenPept'); 
    }
    exit(0);
}
$seq  = $seqio = undef;

eval { 
    ok defined($gb = new Bio::DB::SwissProt('-verbose'=>$verbose)); 
    ok(defined($seq = $gb->get_Seq_by_id('YNB3_YEAST')));
    ok( $seq->length, 125);
    ok($seq->division, 'YEAST');
    ok(defined($seq = $gb->get_Seq_by_acc('P43780')));
    ok( $seq->length, 103); 
    ok( defined( $seq = $gb->get_Seq_by_acc('O39869')));
    ok( $seq->length, 56);
    ok($seq->primary_id, 'O39869');
    ok($seq->division, 'UNK');
 
    # test for bug #958
    $seq = $gb->get_Seq_by_id('P18584');
    ok( defined $seq );
    ok( $seq->length, 497);
    ok( $seq->primary_id, 'DEGP');
    ok( $seq->display_id, 'DEGP_CHLTR');
    ok( $seq->division, 'CHLTR');

    ok( defined($gb = new Bio::DB::SwissProt('-verbose'=>$verbose, 
					     '-retrievaltype' => 'tempfile')));
    ok(defined($seqio = $gb->get_Stream_by_id(['KPY1_ECOLI', 'KPY1_HUMAN'])));
    undef $gb; # testing to see if we can remove gb
    ok( defined($seq = $seqio->next_seq()));
    ok( $seq->length, 470);
    ok( defined($seq = $seqio->next_seq()));
    ok( $seq->length, 530);

};

if ($@) {
    print STDERR "Warning: Couldn't connect to SwissProt with Bio::DB::Swiss.pm!\n$@";
    foreach ( $Test::ntest..$NUMTESTS) { 
	skip(1,1,'could not connect to swissprot');
    }
    exit(0);
}

$seq = undef;

# test the temporary file creation and fasta
eval {
    ok defined ( $gb = new Bio::DB::GenBank('-verbose' =>$verbose,
					    '-format' => 'fasta',
					    '-retrievaltype' => 'tempfile') );
    ok( defined ($seq = $gb->get_Seq_by_id('MUSIGHBA1')));
    ok($seq->length, 408); 
    $seq = $gb->get_Seq_by_acc('AF303112');
    ok( defined $seq);
    ok( $seq->length, 1611);
    # batch mode requires genbank format
    $gb->request_format("genbank");
    ok(defined($seqio = $gb->get_Stream_by_batch([ qw(J00522 AF303112 
							2981014)])));
    ok( $seqio->next_seq->length, 408);
    undef $gb;  # test the case where the db is gone, 
                # but a temp file should remain until seqio goes away. 

    ok($seqio->next_seq->length, 1611);
    ok($seqio->next_seq->length, 1156);
    
};

if ($@) {
    warn "Warning: Couldn't connect to Genbank with Bio::DB::GenBank.pm!\n";
    foreach ( $Test::ntest..$NUMTESTS ) { 
	skip(1,1,'could not connect to Genbank'); 
    }
}
$seq = $seqio = undef;
