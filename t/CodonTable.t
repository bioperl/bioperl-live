#-*-Perl-*-
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
	use vars qw($loaded $DEBUG); }
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

## total number of tests that will be run. 


sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

# create a table object by giving an ID

$myCodonTable = Bio::Tools::CodonTable -> new ( -id => 16);
test 2, defined $myCodonTable && ref($myCodonTable) =~ /Bio::Tools::CodonTable/;

# defaults to ID 1 "Standard"
$myCodonTable = Bio::Tools::CodonTable->new();
test 3, ( $myCodonTable->id() == 1 );


# change codon table
$myCodonTable->id(10);
test 4, defined $myCodonTable;

test 5, ( $myCodonTable->name() eq  'Euplotid Nuclear' );

# translate codons
$myCodonTable->id(1);

@ii  = qw(ACT acu ATN gt ytr sar);
@res = qw(T   T   X   V  L   Z  );
$test = 1;
for $i (0..$#ii) {
    if ($res[$i] ne $myCodonTable->translate($ii[$i]) ) {
	$test = 0; 
	print $ii[$i], ": |", $res[$i], "| ne |", $myCodonTable->translate($ii[$i]), "|\n" if( $DEBUG);
	last ;
    }
}
test 6, $test;

test 7,  ($myCodonTable->translate('ag') eq ''
	  and $myCodonTable->translate('jj') eq ''
	  and $myCodonTable->translate('jjg') eq 'X' 
	  and $myCodonTable->translate('g') eq '');

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
	  $myCodonTable->translate($ii[$i]),  "| @ $i\n" if( $DEBUG);
	last ;
    }
}
test 8, ($test);

# reverse translate amino acids 

@empty = ();
test 9, ($myCodonTable->revtranslate('J') eq @empty);


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
		 print $ii[$i], ': ', $codonres[$j], " ne ", $res[$i][$j], "\n" if( $DEBUG);
		 last TESTING;
	     }
	 }
     }
 }

test 10, ($testing);

#  boolean tests
$myCodonTable->id(1);

test 11, ($myCodonTable->is_start_codon('ATG') and 
	  ! ($myCodonTable->is_start_codon('GGH')) and
	  ($myCodonTable->is_start_codon('HTG')) and
	  ! ($myCodonTable->is_start_codon('CCC')) );

test 12, ($myCodonTable->is_ter_codon('UAG') and  
	  $myCodonTable->is_ter_codon('TaG') and
	  $myCodonTable->is_ter_codon('TaR') and
	  $myCodonTable->is_ter_codon('tRa') and
	  ! ($myCodonTable->is_ter_codon('ttA')));

test 13, ($myCodonTable->is_unknown_codon('jAG') and
	  $myCodonTable->is_unknown_codon('jg') and
	  ! ($myCodonTable->is_unknown_codon('UAG')) );
