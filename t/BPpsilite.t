## Bioperl Test Harness Script for Modules

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
BEGIN { $| = 1; print "1..12\n";
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Tools::BPlite::Iteration;
use Bio::Tools::BPpsilite;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

open FH, "t/psiblastreport.out";
my $report = Bio::Tools::BPpsilite->new(-fh=>\*FH);
test 2, $report;
test 3, $report->query =~ /DICDI/, " query not found";
test 4, $report->database =~ /swissprot/, " database name not found";

my $total_iterations = $report->number_of_iterations;
test 5, $total_iterations == 2, " wrong total iteration number";

my $last_iteration = $report->round($total_iterations);
my $oldhitarray_ref = $last_iteration->oldhits;

# Process initial newly identified hit only
my ($sbjct, $id, $new_hsp, $is_old, @is_old);
 HIT: while($sbjct = $last_iteration->nextSbjct) {
	$id = $sbjct->name;
	$is_old =  grep  /\Q$id\E/, @$oldhitarray_ref;
	if ($is_old ){
		 next HIT;
	}
 	$new_hsp = $sbjct->nextHSP;
	test 6, $new_hsp->score  == 1097, " HSP score not found";
	last HIT;
 }
close FH;

# Verify parsing of PHI-PSI Blast reports
open FH, "t/phipsi.out";
my $report2 = Bio::Tools::BPpsilite->new(-fh=>\*FH);

test 7, $report2;
test 8, $report2->pattern eq "P-E-E-Q", " wrong phi pattern";
test 9, $report2->query_pattern_location->[1] == 120, " wrong phi pattern location";

$total_iterations = $report2->number_of_iterations;
test 10, $total_iterations == 2, " wrong total iteration number in phiblast report";

my $last_iteration2 = $report2->round($total_iterations);
my $sbjct2 = $last_iteration2->nextSbjct;
test 11, $last_iteration2->newhits->[1] =~ /ARATH/, " Hit not found in phiblast report";
my $hsp2 = $sbjct2->nextHSP;
test 12, $hsp2->subject->end == 343, " HSP start not found in phiblast report";

close FH;





