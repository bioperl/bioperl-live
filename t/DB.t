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
BEGIN { $| = 1; print "1..5\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::DB::GenBank;
use Bio::DB::GenPept;
use Bio::DB::SwissProt;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $seq;

eval { 
    my $gb = new Bio::DB::GenBank; 
    $seq = $gb->get_Seq_by_id('MUSIGHBA1'); 
};

if ($@) {
    warn "Warning: Couldn't connect to Genbank with Bio::DB::GenBank.pm!\n{Probably no network access.\n Skipping Test\n";
    warn $@;
    print "ok 2\n";
} else {
    if( $seq->length == 408 ) {
	print "ok 2\n";
    } else {
	print "not ok 2\n";
    }
}
$seq = undef;

eval { 
    my $gb = new Bio::DB::GenPept; 
    $seq = $gb->get_Seq_by_id('195055'); 
};

if ($@) {
    warn "Warning: Couldn't connect to Genbank with Bio::DB::GenPept.pm!\nProbably no network access\n";
    print "ok 3\n";
} else {
    if( $seq->length == 136 ) {
	print "ok 3\n";
    } else {
	print "not ok 3\n";
    }
}

$seq = undef;

eval { 
    my $gb = new Bio::DB::SwissProt; 
    $seq = $gb->get_Seq_by_acc('P43780'); 
};

if ($@) {
    warn ($@);
    warn "Warning: Couldn't connect to Genbank with Bio::DB::Swiss.pm!\nProbably no network access\n";
    print "ok 4\n";
} else {
    if( $seq->length == 103 ) {
	print "ok 4\n";
    } else {
	print "not ok 4\n";
    }
}

print "ok 5\n";




