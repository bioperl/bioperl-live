#-*-perl-*-
# $Id$
use strict;
use warnings;
use Bio::Root::Test;
use Bio::PrimarySeq;

use lib '.';

test_begin(-tests => 91);

use_ok( 'Bio::Restriction::IO' );
use_ok( 'Bio::Restriction::Analysis' );
use_ok('Bio::Restriction::EnzymeCollection');
use_ok('Bio::Restriction::Enzyme');

# recog sites (not nec. cut sites!) in lc

my $seq = new Bio::PrimarySeq(
     -seq        => 'gtcGaagcttAGCAAACGGTTTCTACgacgttatcgtcATTCGGGgcaagcgTCGGCGATTCGGACGTGcacctgcAAAtGCGCGGCgTTAgcgaggtgGCGAgacttttatgtcCCCCTgaagcggttattggTTATATGGTGTTCGTgaccgaTCTAATCCATATTTATTTTTGGCAGTGCtgggtgTTACgacTCGCGA',
     -primary_id => 'test',
     -molecule   => 'dna'
);

# the test enzymes [rebase characterization]: 
# nonambig intrasite cutter: HindIII  [ A^AGCTT        ]
#    ambig intrasite cutter: AasI     [ GACNNNN^NNGTC  ]
# nonambig extrasite cutter: AarI     [ CACCTGC(4/8)   ]
#    ambig extrasite cutter: BceSI    [ SSAAGCG(27/27) ]
#    ambig center    cutter: AjuI     [ (7/12)GAANNNNNNNTTGG(11/6) ]
#    multi extrasite cutter: TaqII    [ GACCGA(11/9),CACCCA(11/9) ]

# the test sequence *cut* (not site) map (recog sites in lc)

#+    AasI(circ)
#+         HindIII                     AasI                      
# 1   CTCGaagcttAGCAAACGGTTTCTACgacgttatcgtcATTCGGGgcaagcgTCGGCGAT  60
#-      IIIdniH                  IsaA      


#+                       BceSI                             AjuI 
#+                       AarI                          AasI
# 61  TCGGACGTGcacctgcAAATGCGCGGCgTTAgcgaggtgGCGAgacttttatgtcCCCCT 120
#-                        IraA                    IsaA            
#-                  ISecB                         IujA      

#+                                                         IIqaT         
#+                             AjuI                TaqII             
# 121 gaagcggttattggTTATATGGTGTTCGTgaccgaTCTAATCCATATTTATTTTTGGCAG 180
#-                    IujA                  IIqaT                 
#-                                                  TaqII          

#+                          
# 181 TGCtgggtgTTACgacTCGCGA 202 
#-              (cric)IsaA   

# so we have +/- cut sites (tpi's, not nt's), in positive strand 
# coordinates:
# HindIII : (5, 9)
# AasI    : (33, 31), (110,108)
# AasI    : (200,202=0) when circularized
# BceSI   : (79, 79)
# AarI    : (80, 84)
# AjuI    : (145, 140) / (113, 108)
# TaqII   : (166, 164) / (174, 172)

ok( my $rebase_io = Bio::Restriction::IO->new(
     -file   => test_input_file('withrefm.906'),
     -format => 'withrefm',
    ), 'read withrefm file');

ok( my $rebase_cln = $rebase_io->read, 'parse withrefm file');

# examples
# ambiguous, nonamibiguous X intrasite, extrasite
ok( my $ninz = $rebase_cln->get_enzyme('HindIII'), 'HindIII: nonambiguous intrasite cutter');
ok( my $nexz = $rebase_cln->get_enzyme('AarI'), 'AarI: nonambiguous extrasite cutter');
ok( my $ainz = $rebase_cln->get_enzyme('AasI'), 'AasI: ambiguous intrasite cutter' );
ok( my $aexz = $rebase_cln->get_enzyme('BceSI'), 'BceSI: ambiguous extrasite cutter' );
# central recognition site: (s/t)[site](m/n)
ok( my $cenz = $rebase_cln->get_enzyme('AjuI'), 'AjuI: cutter with central recog site');
# multisite extrasite:
ok (my $menz = $rebase_cln->get_enzyme('TaqII'), 'TaqII: multi-extrasite cutter');

ok (my $examples = Bio::Restriction::EnzymeCollection->new(
	-enzymes=>[$ninz,$ainz, $ainz, $aexz, $cenz, $menz]
    ) );


# build pretend analysis object to test internals

my $an = {};
bless($an, 'Bio::Restriction::Analysis');
$an->seq($seq);

my ($plus_sites, $minus_sites);
# intrasite cutters

$plus_sites = $an->_make_cuts( $seq->seq, $ninz );
$minus_sites = $an->_make_cuts( $seq->seq, $ninz,'COMP' );
is_deeply( $plus_sites, [5], 'HindIII plus');
is_deeply( $minus_sites, [9], 'HindIII minus');

$plus_sites = $an->_make_cuts(  $seq->seq, $ainz );
$minus_sites = $an->_make_cuts( $seq->seq, $ainz,'COMP' );
is_deeply( $plus_sites, [33, 110], 'AasI plus');
is_deeply( $minus_sites, [31, 108], 'AasI  minus');

# extrasite cutters

$plus_sites = $an->_make_cuts( $seq->seq, $nexz );
$minus_sites = $an->_make_cuts( $seq->seq, $nexz, 'COMP');
is_deeply( $plus_sites, [80], 'AarI plus');
is_deeply( $minus_sites, [84], 'AarI  minus');

$plus_sites = $an->_make_cuts( $seq->seq, $aexz );
$minus_sites = $an->_make_cuts( $seq->seq, $aexz, 'COMP');
is_deeply( $plus_sites, [79], 'BceSI plus');
is_deeply( $minus_sites, [79], 'BceSI  minus');

# central site cutter
$plus_sites = $an->_make_cuts( $seq->seq, $cenz );
$minus_sites = $an->_make_cuts( $seq->seq, $cenz, 'COMP');

is_deeply( $plus_sites, [145, 113], 'AjuI plus');
is_deeply( $minus_sites, [140, 108], 'AjuI minus');

# multisite extrasite cutter
$plus_sites =  $an->_make_cuts( $seq->seq, $menz );
$minus_sites =  $an->_make_cuts( $seq->seq, $menz, 'COMP' );

is_deeply( $plus_sites, [166, 174], 'TaqII plus');
is_deeply( $minus_sites, [164, 172], 'TaqII minus');

# real Analysis object
# start restriction analysis
ok( my $analysis = Bio::Restriction::Analysis->new(
     -seq     => $seq,
     -enzymes => $rebase_cln
    ), "build real B:R::Analysis object");

# retrieve fragment map
my @fm = $analysis->fragment_maps($examples);
is( @fm, 13, '13 fragments');
# circularize
ok( $seq->is_circular(1), 'circularize');
ok( $analysis->cut, 'recut');
@fm = $analysis->fragment_maps($examples);
is_deeply( [$analysis->positions('AasI')], [33, 110, 200], 'circ: AasI
site at origin' );
is( @fm, 13, 'circ: still 13 fragments (cut site at origin)');


# Emmanuel's tests / bug3015

use_ok('Bio::Restriction::IO');
use_ok('Bio::Restriction::Analysis');

#ATGAGCGCTcacgtcACTAG^TTCTAGAGcctcagcAATTC^CGATCccgctcGATTAATGC^TccgcAGCAGCGATATCGAG^CATGGTCATGAgaatgcGGC^ATCGATCGGCATTATATcacgtcAATCGCGTCGCTGCATGCTAGCG
$seq = new Bio::PrimarySeq(
     -seq        => 'ATGAGCGCTcacgtcACTAGTTCTAGAGcctcagcAATTCCGATCccgctcGATTAATGCTccgcAGCAGCGATATCGAGCATGGTCATGAgaatgcGGCATCGATCGGCATTATATcacgtcAATCGCGTCGCTGCATGCTAGCG',
     -primary_id => 'test1',
     -molecule   => 'dna'
);

ok( $rebase_io = Bio::Restriction::IO->new(
     -file   => test_input_file('withrefm.906'),
     -format => 'withrefm',
    ), 'read withrefm file');

ok( $rebase_cln = $rebase_io->read, 'parse withrefm file');

my @enzs = qw/AbeI AccBSI AciI Asp26HI BmgBI/;

# Get the enzymes with back cut site
# 'AbeI'    CCTCAGC(-5/-2) #1
# 'AccBSI'  CCGCTC(-3/-3)  #1
# 'AciI'    CCGC(-3/-1)    #2
# 'Asp26HI' GAATGC(1/-1)   #1
# 'BmgBI'   CACGTC(-3/-3)  #2

ok(my $collection = Bio::Restriction::EnzymeCollection->new(-empty => 1), "Collection initiated");

foreach my $e (@enzs){

    ok(my $enz = $rebase_cln->get_enzyme($e), "$e: found ok into collection");
    $collection->enzymes($enz);
}

$an = { };
bless($an, 'Bio::Restriction::Analysis');
$an->seq($seq);
$an->enzymes($collection);

#Test all types of stuff from the enzyme
my $data = {
   'AbeI'    => { 'p' => [23],
                  'm' => [26],
                  'pos' => 2,
                  'f' => 3,
                  'o' => "5'",
                  's' => 'CC^TCAGC',
                  'r' => 'GCTGAGG',
                  'c' => '-5',
                  'rc' => '-2' },

   'AccBSI'  => { 'p' => [42],
                  'm' => [42],
                  'pos' => 1,
                  'f' => 2,
                  'o' => "blunt",
                  's' => 'CCG^CTC',
                  'r' => 'GAGCGG',
                  'c' => '-3',
                  'rc' => '-3' },

   'AciI'    => { 'p' => [42,58,100],
                  'm' => [44,60,102],
                  'pos' => 6,
                  'f' => 7,
                  'o' => "5'",
                  's' => 'C^CGC',
                  'r' => 'GCGG',
                  'c' => '-3',
                  'rc' => '-1' },

   'Asp26HI' => { 'p' => [98],
                  'm' => [90],
                  'pos' => 1,
                  'f' => 2,
                  'o' => "3'",
                  's' => 'GAATGC',
                  'r' => 'GCATTC',
                  'c' => 7,
                  'rc' => '-1' },

   'BmgBI'   => { 'p' => [6, 114],
                  'm' => [6, 114],
                  'pos' => 2,
                  'f' => 3,
                  'o' => "blunt",
                  's' => 'CAC^GTC',
                  'r' => 'GACGTG',
                  'c' => '-3',
                  'rc' => '-3' },
};

foreach my $e (@enzs){

   my $d = $data->{$e};
   my $z = $rebase_cln->get_enzyme($e);
   my $minus_sites = $an->_make_cuts($seq->seq, $z, 'COMP');
   my $plus_sites  = $an->_make_cuts($seq->seq, $z);

   is_deeply($plus_sites,  $d->{p}, "$e plus");
   is_deeply($minus_sites, $d->{m}, "$e minus");

   ok(scalar($an->fragments($z)) eq $d->{f}, "$e fragment");
   ok(scalar($an->positions($e)) eq $d->{pos}, "$e positions");

   is($z->overhang(), $d->{'o'}, "$e Overhang");
   is($z->name(), "$e", "$e name");
   is($z->site(), $d->{s}, "$e site");
   is($z->revcom_site(), $d->{r}, "$e revcom_site");
   is($z->cut(), $d->{c}, "$e cut");
   is($z->complementary_cut(), $d->{rc}, "$e complementary_cut");
}

