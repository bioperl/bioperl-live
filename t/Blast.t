## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
# Test script for Bio::Tools::Blast.pm
# Steve A. Chervitz, sac@genome.stanford.edu
# Fairly rudimentary. Only tests parsing code.
# Strategy for this test was borrowed from L. Stein's BoulderIO package.

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

BEGIN {	
    $| = 1; print "1..29\n"; 
    $testout = "blast.t.out";  # output from this script.
    $expectedout = 't/expected.blast.out';
    unlink $testout;
    $^W = 0; 
}
END {
#   print "not ok 1\n" unless $loaded;
#   unlink $testout;  # commented out since you may want to check it...
}

use lib '../';
use Bio::Tools::Blast;

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}
my($s,@s);

open (OUT,">$testout");

test 1, $blast = Bio::Tools::Blast->new(-file   =>'t/blast.report',
					-signif => 1e-5,
					-parse  => 1,
					-stats  => 1,
					-check_all_hits => 1,
					);
test 2, $blast->display();
test 3, $blast->is_signif;
test 4, $blast->signif eq '1.0e-05', "Signif: ".$blast->signif;
test 5, $blast->num_hits == 4;
test 6, $blast->length == 504;
test 7, $blast->program eq 'TBLASTN';
test 8, $blast->query eq 'gi|1401126';
test 9, $blast->hit->name eq 'gb|U49928|HSU49928';
#print STDERR "Hit is ",$blast->hit->name,"\n";
test 10, $blast->hit->length == 3096;

@hits  = $blast->hits;

test 11, $hits[0]->expect eq '0.0';
test 12, $hits[1]->expect eq '4e-07';
test 13, $hits[2]->expect eq '1e-05';
test 14, $hits[1]->frac_identical eq '0.25';
test 15, $hits[1]->hsp->frac_conserved eq '0.43';
test 16, $hits[1]->hsp->score == 137;
test 17, $hits[1]->hsp->bits eq '57.8';

# Sequence index testing.
test 18, @inds = $hits[1]->hsp->seq_inds('query', 'iden', 1);
test 19, $inds[0] eq '66-68';

# Output testing.
test 20, print OUT $blast->table_labels;
test 21, print OUT $blast->table;
print OUT "\n\n";
test 22, print OUT $blast->table_labels_tiled;
test 23, print OUT $blast->table_tiled;
print OUT "\n\n";
close OUT;

test 24, -s $blast->file;
test 25, $cfile = $blast->compress_file;
test 26, -s $cfile and -B $cfile, "Can't compress Blast file";
test 27, $ufile = $blast->uncompress_file;
test 28, -s $ufile and -T $ufile, "Can't uncompress Blast file";

print "checking expected output...\n";

test 29, system('diff', $testout, $expectedout) == 0, 
    "diff $testout $expectedout";







