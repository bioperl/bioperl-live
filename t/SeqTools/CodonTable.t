# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use Bio::Root::Test;

    test_begin(-tests => 94);

    use_ok('Bio::Tools::CodonTable');
    use_ok('Bio::CodonUsage::IO');
}

# create a table object by giving an ID
my $DEBUG = test_debug();
my $myCodonTable = Bio::Tools::CodonTable -> new ( -id => 16);
ok defined $myCodonTable;
isa_ok $myCodonTable, 'Bio::Tools::CodonTable';

# Access to ID table 0 through constructor
$myCodonTable = Bio::Tools::CodonTable->new( -id => 0);
is $myCodonTable->id(), 0;

# defaults to ID 1 "Standard"
$myCodonTable = Bio::Tools::CodonTable->new();
is $myCodonTable->id(), 1;
is $myCodonTable->name(), "Standard";


# Test invalid IDs warn and return table id=1
{
    # ID table 7 is invalid because it has been removed
    foreach my $invalid_id (99, -2, 7) {
        my $stderr = '';
        # capture stderr output
        local *STDERR;
        open STDERR, '>', \$stderr;
        $myCodonTable = Bio::Tools::CodonTable->new(-id => $invalid_id);
        like $stderr, qr/Not a valid codon table ID/;
        is $myCodonTable->id, 1;
    }
}


# change codon table
$myCodonTable->id(10);
is $myCodonTable->id, 10;
is $myCodonTable->name(), 'Euplotid Nuclear';

# enumerate tables as object method
my $table = $myCodonTable->tables();
cmp_ok (keys %{$table}, '>=', 26); # currently 25 known tables + "strict"
is $table->{11}, 'Bacterial, Archaeal and Plant Plastid';

# enumerate tables as class method
$table = Bio::Tools::CodonTable->tables;
cmp_ok (values %{$table}, '>=', 26); # currently 25 known tables + "strict"
is $table->{23}, 'Thraustochytrium Mitochondrial';

# translate codons
$myCodonTable->id(1);

eval {
    $myCodonTable->translate();
};
ok ($@ =~ /EX/) ;

# Automatically completing translation of incomplete codons is no longer default
# behavior b/c of inconsistent behavior compared with Bio::PrimarySeq::translate
# and unexpected side effects (e.g. what if the last few bases isn't supposed to
# be translated). To re-establish this, pass a second argument to the method.

is $myCodonTable->translate(''), '';

my @ii  = qw(ACT acu ATN gt ytr sar);
my @res = qw(T   T   X   V  L   Z  );
my $test = 1;
for my $i (0..$#ii) {
    if ($res[$i] ne $myCodonTable->translate($ii[$i], 1) ) {
        $test = 0;
        print $ii[$i], ": |", $res[$i], "| ne |",
        $myCodonTable->translate($ii[$i], 1), "|\n" if( $DEBUG);
        last ;
    }
}
ok ($test);
is $myCodonTable->translate('ag'), '';
is $myCodonTable->translate('ag',1), '';

is $myCodonTable->translate('jj'), '';
is $myCodonTable->translate('jj',1), '';

is $myCodonTable->translate('jjg'), 'X';
is $myCodonTable->translate('jjg',1), 'X';

is $myCodonTable->translate('gt'), '';
is $myCodonTable->translate('gt',1), 'V';

is $myCodonTable->translate('g'), '';
is $myCodonTable->translate('g',1), '';

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

is $myCodonTable->revtranslate('U'), 0;
is $myCodonTable->revtranslate('O'), 0;
is $myCodonTable->revtranslate('J'), 9;
is $myCodonTable->revtranslate('I'), 3;
my @RNA_codons = $myCodonTable->revtranslate('M', 'RNA');
is $RNA_codons[0], 'aug'; # test RNA output

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

# boolean tests
$myCodonTable->id(1); # Standard table

ok( $myCodonTable->is_start_codon('ATG'), 'is_start_codon, ATG');
ok( ! $myCodonTable->is_start_codon('GGH'), 'is_start_codon, GGH');
ok( $myCodonTable->is_start_codon('HTG'), 'is_start_codon, HTG');
ok( ! $myCodonTable->is_start_codon('CCC'), 'is_start_codon, CCC');

ok( $myCodonTable->is_ter_codon('UAG'), 'is_ter_codon, U should map to T, UAG');
ok( $myCodonTable->is_ter_codon('TaG'), 'is_ter_codon,TaG');
ok( $myCodonTable->is_ter_codon('TaR'), 'is_ter_codon,TaR');
ok( $myCodonTable->is_ter_codon('tRa'), 'is_ter_codon,tRa');
ok( ! $myCodonTable->is_ter_codon('ttA'), 'is_ter_codon,ttA');

# Ambiguous codons should fail
ok ! $myCodonTable->is_ter_codon('NNN'), 'is_ter_codon, ambiguous codons should fail, NNN';
ok ! $myCodonTable->is_ter_codon('TAN'), 'is_ter_codon, ambiguous codons should fail, TAN';
ok ! $myCodonTable->is_ter_codon('CC'), 'is_ter_codon, incomplete codons should fail, CC';

ok ! $myCodonTable->is_start_codon('NNN'), 'is_start_codon, ambiguous codons should fail, NNN';
ok ! $myCodonTable->is_start_codon('NTG'), 'is_start_codon, ambiguous codons should fail, NTG';
ok ! $myCodonTable->is_start_codon('N'), 'is_start_codon, incomplete codons should fail, NN';

ok $myCodonTable->is_unknown_codon('jAG');
ok $myCodonTable->is_unknown_codon('jg');
ok ! $myCodonTable->is_unknown_codon('UAG');

is $myCodonTable->translate_strict('ATG'), 'M';

#
# adding a custom codon table
#

my @custom_table =
    ( 'test1',
      'FFLLSSSSYY**CC*WLLLL**PPHHQQR*RRIIIMT*TT*NKKSSRRV*VVAA*ADDEE*GGG'
    );

ok my $custct = $myCodonTable->add_table(@custom_table);
is $custct, 33;
is $myCodonTable->translate('atgaaraayacmacracwacka'), 'MKNTTTT';
ok $myCodonTable->id($custct);
is $myCodonTable->translate('atgaaraayacmacracwacka'), 'MKXXTTT';

# test doing this via Bio::PrimarySeq object

use Bio::PrimarySeq;
ok $seq = Bio::PrimarySeq->new(-seq=>'atgaaraayacmacracwacka', -alphabet=>'dna');
is $seq->translate()->seq, 'MKNTTTT';
is $seq->translate(undef, undef, undef, undef, undef, undef, $myCodonTable)->seq, 'MKXXTTT';

# test gapped translated

ok $seq = Bio::PrimarySeq->new(-seq      => 'atg---aar------aay',
                               -alphabet => 'dna');
is $seq->translate->seq, 'M-K--N';

ok $seq = Bio::PrimarySeq->new(-seq =>'ASDFGHKL');
is $myCodonTable->reverse_translate_all($seq), 'GCBWSNGAYTTYGGVCAYAARYTN';
ok $seq = Bio::PrimarySeq->new(-seq => 'ASXFHKL');
is $myCodonTable->reverse_translate_all($seq), 'GCBWSNNNNTTYCAYAARYTN';

#
# test reverse_translate_best(), requires a Bio::CodonUsage::Table object
#

ok $seq = Bio::PrimarySeq->new(-seq =>'ACDEFGHIKLMNPQRSTVWYX');
ok my $io = Bio::CodonUsage::IO->new(-file => test_input_file('MmCT'));
ok my $cut = $io->next_data();
is $myCodonTable->reverse_translate_best($seq,$cut), 'GCCTGCGACGAGTTCGGCCACATCAAGCTGATGAACCCCCAGCGCTCCACCGTGTGGTACNNN';
is $myCodonTable->reverse_translate_all($seq, $cut, 15), 'GCNTGYGAYGARTTYGGVCAYATYAARCTSATGAAYCCNCARMGVWSYACHGTSTGGTAYNNN';

#
# test 'Strict' table, requires a Bio::CodonUsage::Table object
#

$myCodonTable = Bio::Tools::CodonTable->new(); # Default Standard table

#  boolean tests
ok( $myCodonTable->is_start_codon('ATG'));
ok( ! $myCodonTable->is_start_codon('GTG'));
ok( $myCodonTable->is_start_codon('TTG'));
ok( $myCodonTable->is_start_codon('CTG'));
ok( ! $myCodonTable->is_start_codon('CCC'));

$myCodonTable->id(0); # Special 'Strict' table (ATG-only start)
ok( $myCodonTable->is_start_codon('ATG'));
ok( ! $myCodonTable->is_start_codon('GTG'));
ok( ! $myCodonTable->is_start_codon('TTG'));
ok( ! $myCodonTable->is_start_codon('CTG'));
ok( ! $myCodonTable->is_start_codon('CCC'));

# Pterobranchia Mitochondrial codon table
$myCodonTable->id(24);
ok( $myCodonTable->is_start_codon('GTG'));
ok( $myCodonTable->is_start_codon('CTG'));
is( $myCodonTable->translate_strict('TGA'), 'W');

# Candidate Division SR1 and Gracilibacteria codon table
$myCodonTable->id(25);
ok( $myCodonTable->is_start_codon('GTG'));
ok( !$myCodonTable->is_start_codon('CTG'));
is($myCodonTable->translate_strict('TGA'), 'G');

# The name for this is pretty long and spans mutliple lines.  Confirm
# it was parsed correctly from the ASN.1 file.
$myCodonTable->id(4);
is $myCodonTable->name(), "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma";
