# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { 
    eval { require 'IO/String.pm' };
    if( $@ ) {
	print STDERR "IO::String not loaded. This means DB test cannot be executed. Skipping\n";
	print "1..1\n";
	print "ok 1\n";
	exit(0);
    } 
    
    $| = 1; print "1..25\n"; 
    use vars qw($loaded $testnum); 
}

END {print "not ok 1\n" unless $loaded;}

use Bio::DB::GenBank;
use Bio::DB::GenPept;
use Bio::DB::SwissProt;
use strict;
$loaded = 1;
my $testnum;
my $verbose = 0;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

{ 
    $testnum = 2;
    sub test ($;$) {
	my($true,$msg) = @_;	
	$msg = '' if !defined $msg;
	print($true ? "ok $testnum\n" : "not ok $testnum $msg\n");
	$testnum++;
    }
}

my ($gb,$seq,$seqio);
# get a single seq
eval { 
    test defined ( $gb = new Bio::DB::GenBank(-verbose=>$verbose) );     
    test defined ($seq = $gb->get_Seq_by_id('MUSIGHBA1'))
	&& $seq->length == 408; 
    test  defined ($seq = $gb->get_Seq_by_acc('AF303112')) 
	&& $seq->length == 1611; 
    test defined($seqio = $gb->get_Stream_by_batch([ qw(J00522 AF303112 
							   2981014)]));
    test $seqio->next_seq->length == 408;
    test $seqio->next_seq->length == 1611;
    test $seqio->next_seq->length == 1156;
};

if ($@) {
    warn "Warning: Couldn't connect to Genbank with Bio::DB::GenBank.pm!\nProbably no network access.\n Skipping Test\n";
    warn $@;
    while ( $testnum <= 8 ) { test 0 }
}
$seq = $seqio = undef;

eval { 
    test defined($gb = new Bio::DB::GenPept(-verbose=>$verbose)); 
    test defined($seq = $gb->get_Seq_by_id('195055'))
	&& $seq->length == 136; 
    test defined($seq = $gb->get_Seq_by_acc('AAC06201'))
	&& $seq->length == 353;
    test defined($seqio = $gb->get_Stream_by_batch([ qw(AAC06201 195055)]));
    test $seqio->next_seq->length() == 353;
    test $seqio->next_seq->length == 136;
};

if ($@) {
    warn "Warning: Couldn't connect to Genbank with Bio::DB::GenPept.pm!\nProbably no network access\n";
    warn $@;
    while( $testnum <= 14 ) { test 0 }    
}
$seq  = $seqio = undef;

eval { 
    test defined($gb = new Bio::DB::SwissProt(-verbose=>$verbose)); 
    test defined($seq = $gb->get_Seq_by_acc('P43780')) 
	&& $seq->length == 103; 
    $gb = new Bio::DB::SwissProt(-verbose=>$verbose, 
				 -retrievaltype => 'tempfile');
    test defined($seqio = $gb->get_Stream_by_id('KPY1_ECOLI'));
    undef $gb; # testing to see if we can remove gb
    test defined($seq = $seqio->next_seq()) && $seq->length == 470;
};

if ($@) {
    print STDERR "Warning: Couldn't connect to Genbank with Bio::DB::Swiss.pm!\nProbably no network access\n" . $@;

    while( $testnum <= 18) { test 0;}

}
$seq = undef;

# test the temporary file creation and fasta
eval {
    test defined ( $gb = new Bio::DB::GenBank(-verbose=>$verbose,
					      -format => 'fasta',
					      -retrievaltype => 'tempfile') );
    test defined ($seq = $gb->get_Seq_by_id('MUSIGHBA1'))
	&& $seq->length == 408; 
    test  defined ($seq = $gb->get_Seq_by_acc('AF303112')) 
	&& $seq->length == 1611; 
    test defined($seqio = $gb->get_Stream_by_batch([ qw(J00522 AF303112 
							   2981014)]));
    test $seqio->next_seq->length == 408;
    undef $gb;  # test the case where the db is gone, 
                # but a temp file should remain until seqio goes away. 

    test $seqio->next_seq->length == 1611;
    test $seqio->next_seq->length == 1156;
    
};

if ($@) {
    warn "Warning: Couldn't connect to Genbank with Bio::DB::GenBank.pm!\nProbably no network access.\n Skipping Test\n";
    warn $@;
    while ( $testnum <= 25 ) { test 0 }
}
$seq = $seqio = undef;
