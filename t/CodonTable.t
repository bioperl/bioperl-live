## Bioperl Test Harness Script for Modules
## $Id$

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
BEGIN { $| = 1; print "1..10\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

#use lib '../';
use Bio::Tools::CodonTable;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


# create a table object by giving an ID

$myCodonTable = Bio::Tools::CodonTable -> new ( -id => 16);
print "ok 2\n";  


#($start, $stop) = $range->intersection($range2);
#if( $start == 15 && $stop == 20 ) {
#   print "ok 5\n";
#} else {
#   print "not ok 5\n";
#}      

# defaults to ID 1 "Standard"
$myCodonTable = Bio::Tools::CodonTable->new();
if ( $myCodonTable->id() == 1 ) {
   print "ok 3\n";
} else {
   print "not ok 3\n";
}

#check if a codon table ID is valid before setting it:
if ( $myCodonTable->id(55) ) {
   print "not ok 4\n";
}
else {
   print "ok 4\n";
}

# change codon table
$myCodonTable->id(10);
print "ok 5\n";

if( $myCodonTable->name() eq  'Euplotid Nuclear' ) {
   print "ok 6\n";
} else {
   print "not ok 6\n";
}

# translate codons
@ii = ('ACT', 'acu', 'AC', 'ACN');
@res = ('T', 'T', '', 'X' );
foreach $i (@ii) {
    $aa = $myCodonTable->translate($i);
    push @aa, $aa;
}
if (@aa == @res) {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
}


# reverse translate amino acids 
@aa=();
@ii = ('A',  'g',  'ACN', 'Thr', 'sER' );
@res = qw(gct gcc gca gcg ggt ggc gga ggg act acc aca acg tct tcc tca tcg agt agc);
foreach $i (@ii) {
    @codon1 = $myCodonTable->revtranslate($i);
    push (@aa, @codon1);
}

if (@aa == @res) {
    print "ok 8\n";
}
else {
    print "not ok 8\n";
}

#  #boolean tests
  $myCodonTable->id(1);

if ($myCodonTable->is_start_codon('ATG') and 
    !  ($myCodonTable->is_start_codon('GGH')) and
    !  ($myCodonTable->is_start_codon('CCC')) ){
    print "ok 9\n";
}
else {
    print "not ok 9\n";
}

if ($myCodonTable->is_ter_codon('UAG') and  
    $myCodonTable->is_ter_codon('TaG') and
    ! ($myCodonTable->is_ter_codon('ttA'))) {
    print "ok 10\n";
}
else {
    print "not ok 10\n";
}
