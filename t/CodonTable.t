#-*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { 
	use lib 't'; 
    }
    use Test;

    plan tests => 28;
}
use Bio::Tools::CodonTable;
use vars qw($DEBUG);
ok(1);

# create a table object by giving an ID

my $myCodonTable = Bio::Tools::CodonTable -> new ( -id => 16);
ok defined $myCodonTable;
ok $myCodonTable->isa('Bio::Tools::CodonTable');

# defaults to ID 1 "Standard"
$myCodonTable = Bio::Tools::CodonTable->new();
ok $myCodonTable->id(), 1;


# change codon table
$myCodonTable->id(10);
ok $myCodonTable->id, 10;

ok $myCodonTable->name(), 'Euplotid Nuclear';

# translate codons
$myCodonTable->id(1);

my @ii  = qw(ACT acu ATN gt ytr sar);
my @res = qw(T   T   X   V  L   Z  );
my $test = 1;
for my $i (0..$#ii) {
    if ($res[$i] ne $myCodonTable->translate($ii[$i]) ) {
	$test = 0; 
	print $ii[$i], ": |", $res[$i], "| ne |", $myCodonTable->translate($ii[$i]), "|\n" if( $DEBUG);
	last ;
    }
}
ok $test;

ok $myCodonTable->translate('ag'), '';
ok $myCodonTable->translate('jj'), '';
ok $myCodonTable->translate('jjg'), 'X';
ok $myCodonTable->translate('gt'), 'V'; 
ok $myCodonTable->translate('g'), '';

# a more comprehensive test on ambiguous codes
my $seq = <<SEQ;
atgaaraayacmacracwackacyacsacvachacdacbacxagragyatmatwatyathcarcayc
cmccrccwcckccyccsccvcchccdccbccxcgmcgrcgwcgkcgycgscgvcghcgdcgbcgxctmctrct
wctkctyctsctvcthctdctbctxgargaygcmgcrgcwgckgcygcsgcvgchgcdgcbgcxggmggrggw
ggkggyggsggvgghggdggbggxgtmgtrgtwgtkgtygtsgtvgthgtdgtbgtxtartaytcmtcrtcwt
cktcytcstcvtchtcdtcbtcxtgyttrttytramgamggmgrracratrayytaytgytrsaasagsartaa;
SEQ
    $seq =~ s/\s+//g;
@ii = grep { length == 3 } split /(.{3})/, $seq; 
print join (' ', @ii), "\n" if( $DEBUG);
my $prot = <<PROT;
MKNTTTTTTTTTTTRSIIIIQHPPPPPPPPPPPRRRRRRRRRRRLLLLLLLLLLLEDAAAAAAAAAAAGGG
GGGGGGGGVVVVVVVVVVV*YSSSSSSSSSSSCLF*RRRBBBLLLZZZ*
PROT

    $prot =~ s/\s//;
@res = split //, $prot;
print join (' ', @res), "\n" if( $DEBUG );
$test = 1;
for my $i (0..$#ii) {
    if ($res[$i] ne $myCodonTable->translate($ii[$i]) ) {
	$test = 0; 
	print $ii[$i], ": |", $res[$i], "| ne |", 
	  $myCodonTable->translate($ii[$i]),  "| @ $i\n" if( $DEBUG);
	last ;
    }
}
ok $test;

# reverse translate amino acids 

ok $myCodonTable->revtranslate('J'), 0;


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

$test = 1;
 TESTING: {
     for my $i (0..$#ii) {
	 my @codonres = $myCodonTable->revtranslate($ii[$i]);
	 for my $j (0..$#codonres) {
	     if ($codonres[$j] ne $res[$i][$j]) {
		 $test = 0;
		 print $ii[$i], ': ', $codonres[$j], " ne ", 
		 $res[$i][$j], "\n" if( $DEBUG);
		 last TESTING;
	     }
	 }
     }
 }
ok $test;

#  boolean tests
$myCodonTable->id(1);

ok $myCodonTable->is_start_codon('ATG');  
ok $myCodonTable->is_start_codon('GGH'), 0;
ok $myCodonTable->is_start_codon('HTG');
ok $myCodonTable->is_start_codon('CCC'), 0;

ok $myCodonTable->is_ter_codon('UAG');
ok $myCodonTable->is_ter_codon('TaG');
ok $myCodonTable->is_ter_codon('TaR');
ok $myCodonTable->is_ter_codon('tRa');
ok $myCodonTable->is_ter_codon('ttA'), 0;

ok $myCodonTable->is_unknown_codon('jAG');
ok $myCodonTable->is_unknown_codon('jg');
ok $myCodonTable->is_unknown_codon('UAG'), 0;


ok $myCodonTable->translate_strict('ATG'), 'M';
