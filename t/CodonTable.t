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
BEGIN { $| = 1; print "1..13\n"; 
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


# defaults to ID 1 "Standard"
$myCodonTable = Bio::Tools::CodonTable->new();
if ( $myCodonTable->id() == 1 ) {
   print "ok 3\n";
} else {
   print "not ok 3\n";
}


# change codon table
$myCodonTable->id(10);
print "ok 4\n";

if( $myCodonTable->name() eq  'Euplotid Nuclear' ) {
   print "ok 5\n";
} else {
   print "not ok 5\n";
}

# translate codons
$myCodonTable->id(1);

@ii  = qw(ACT acu ATN gt ytr sar);
@res = qw(T   T   X   V  L   Z  );
$test = 1;
for $i (0..$#ii) {
    if ($res[$i] ne $myCodonTable->translate($ii[$i]) ) {
	$test = 0; 
	print $ii[$i], ": |", $res[$i], "| ne |", $myCodonTable->translate($ii[$i]), "|\n";
	last ;
    }
}
if ($test) {
    print "ok 6\n";
}
else {
    print "not ok 6\n";
}


if ($myCodonTable->translate('ag') eq ''
    and $myCodonTable->translate('jj') eq ''
    and $myCodonTable->translate('jjg') eq 'X' 
    and $myCodonTable->translate('g') eq '') {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
}



# a more comprehensive test on ambiguous codes
$seq = "atgaaraayacmacracwackacyacsacvachacdacbacxagragyatmatwatyathcarcayc".
    "cmccrccwcckccyccsccvcchccdccbccxcgmcgrcgwcgkcgycgscgvcghcgdcgbcgxctmctrct".
    "wctkctyctsctvcthctdctbctxgargaygcmgcrgcwgckgcygcsgcvgchgcdgcbgcxggmggrggw".
    "ggkggyggsggvgghggdggbggxgtmgtrgtwgtkgtygtsgtvgthgtdgtbgtxtartaytcmtcrtcwt".
    "cktcytcstcvtchtcdtcbtcxtgyttrttytramgamggmgrracratrayytaytgytrsaasagsartaa";

@ii = grep { length == 3 } split /(.{3})/, $seq; 
#print join (' ', @ii), "\n"; 
$prot = 'MKNTTTTTTTTTTTRSIIIIQHPPPPPPPPPPPRRRRRRRRRRRLLLLLLLLLLLEDAAAAAAAAAAAGGG'.
    'GGGGGGGGVVVVVVVVVVV*YSSSSSSSSSSSCLF*RRRBBBLLLZZZ*';
@res = split //, $prot;
#print join (' ', @res), "\n";
$test = 1;
for $i (0..$#ii) {
    if ($res[$i] ne $myCodonTable->translate($ii[$i]) ) {
	$test = 0; 
	print $ii[$i], ": |", $res[$i], "| ne |", 
	  $myCodonTable->translate($ii[$i]),  "| @ $i\n";
	last ;
    }
}
if ($test) {
    print "ok 8\n";
}
else {
    print "not ok 8\n";
}

# reverse translate amino acids 

@empty = ();
if ($myCodonTable->revtranslate('J') eq @empty) {
    print "ok 9\n";
}
else {
    print "not ok 9\n";
}


@ii = qw(A l ACN Thr sER ter Glx);
@res = (
	[qw(gct gcc gca gcg)],
	[qw(ggc gga ggg act acc aca acg)],
	[qw(tct tcc tca tcg agt agc)],
	[qw(act acc aca acg)],
	[qw(tct tcc tca tcg agt agc)],
	[qw(taa tag tga)],
	[qw(gaa gag caa cag)]
	);

$testing = 1;
 TESTING: {
     for $i (0..$#ii) {
	 @codonres = $myCodonTable->revtranslate($ii[$i]);
	 for $j (0..$#codonres) {
	     if ($codonres[$j] ne $res[$i][$j]) {
		 $testing = 0;
		 print $ii[$i], ': ', $codonres[$j], " ne ", $res[$i][$j], "\n";
		 last TESTING;
	     }
	 }
     }
 }

if ($testing) {
    print "ok 10\n";
}
else {
    print "not ok 10\n";
}

#  boolean tests
$myCodonTable->id(1);

if ($myCodonTable->is_start_codon('ATG') and 
    ! ($myCodonTable->is_start_codon('GGH')) and
    ($myCodonTable->is_start_codon('HTG')) and
    ! ($myCodonTable->is_start_codon('CCC')) ){
    print "ok 11\n";
}
else {
    print "not ok 11\n";
}

if ($myCodonTable->is_ter_codon('UAG') and  
    $myCodonTable->is_ter_codon('TaG') and
    $myCodonTable->is_ter_codon('TaR') and
    $myCodonTable->is_ter_codon('tRa') and
    ! ($myCodonTable->is_ter_codon('ttA'))) {
    print "ok 12\n";
}
else {
    print "not ok 12\n";
}

if ($myCodonTable->is_unknown_codon('jAG') and
    $myCodonTable->is_unknown_codon('jg') and
    ! ($myCodonTable->is_unknown_codon('UAG')) ) {  
    print "ok 13\n";
}
else {
    print "not ok 13\n";
}
