# -*-Perl-*-
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
    plan tests => 109;
}

use Bio::Phenotype::Phenotype;
use Bio::Species;
use Bio::Annotation::Reference;
use Bio::Map::CytoPosition;
use Bio::Phenotype::Correlate;
use Bio::Phenotype::Measure;
use Bio::Annotation::DBLink;


my $obj = Bio::Phenotype::Phenotype->new();

ok( $obj->isa( "Bio::Phenotype::PhenotypeI" ) );
ok( $obj->isa( "Bio::Phenotype::Phenotype" ) );

ok( $obj->name( "r1" ) );
ok( $obj->name(), "r1" );

ok( $obj->description( "This is ..." ) );
ok( $obj->description(), "This is ..." );

my $mouse = Bio::Species->new();
$mouse->classification( qw( musculus Mus ) );
ok( $obj->species( $mouse ) );
ok( $obj->species()->binomial(), "Mus musculus" );

ok( $obj->comment( "putative" ) );
ok( $obj->comment(), "putative" );



ok( $obj->each_gene_symbol(), 0 );

ok( $obj->add_gene_symbols( ( "A", "B" ) ) );
ok( $obj->each_gene_symbol(), 2 );
my @gs = $obj->each_gene_symbol();
ok( $gs[ 0 ], "A" );
ok( $gs[ 1 ], "B" );
ok( $obj->each_gene_symbol(), 2 );

my @gs2 = $obj->remove_gene_symbols();
ok( $gs2[ 0 ], "A" );
ok( $gs2[ 1 ], "B" );

ok( $obj->each_gene_symbol(), 0 );
ok( $obj->remove_gene_symbols(), 0 );



my $v1 = Bio::Variation::VariantI->new();
my $v2 = Bio::Variation::VariantI->new();

$v1->length( "123" );

ok( $obj->each_Variant(), 0 );

ok( $obj->add_Variants( ( $v1, $v2 ) ) );
ok( $obj->each_Variant(), 2 );
my @vs = $obj->each_Variant();
ok( $vs[ 0 ], $v1 );
ok( $vs[ 1 ], $v2 );
ok( $vs[ 0 ]->length(), "123" );
ok( $obj->each_Variant(), 2 );

my @vs2 = $obj->remove_Variants();
ok( $vs2[ 0 ], $v1 );
ok( $vs2[ 1 ], $v2 );

ok( $obj->each_Variant(), 0 );
ok( $obj->remove_Variants(), 0 );




my $r1 = Bio::Annotation::Reference->new();
my $r2 = Bio::Annotation::Reference->new();

$r1->title( "title" );

ok( $obj->each_Reference(), 0 );

ok( $obj->add_References( ( $r1, $r2 ) ) );
ok( $obj->each_Reference(), 2 );
my @rs = $obj->each_Reference();
ok( $rs[ 0 ], $r1 );
ok( $rs[ 1 ], $r2 );
ok( $rs[ 0 ]->title(), "title" );
ok( $obj->each_Reference(), 2 );

my @rs2 = $obj->remove_References();
ok( $rs2[ 0 ], $r1 );
ok( $rs2[ 1 ], $r2 );

ok( $obj->each_Reference(), 0 );
ok( $obj->remove_References(), 0 );




my $c1 = Bio::Map::CytoPosition->new();
my $c2 = Bio::Map::CytoPosition->new();

$c1->chr( "12" );

ok( $obj->each_CytoPosition(), 0 );

ok( $obj->add_CytoPositions( ( $c1, $c2 ) ) );
ok( $obj->each_CytoPosition(), 2 );
my @cs = $obj->each_CytoPosition();
ok( $cs[ 0 ], $c1 );
ok( $cs[ 1 ], $c2 );
ok( $cs[ 0 ]->chr(), "12" );
ok( $obj->each_CytoPosition(), 2 );

my @cs2 = $obj->remove_CytoPositions();
ok( $cs2[ 0 ], $c1 );
ok( $cs2[ 1 ], $c2 );

ok( $obj->each_CytoPosition(), 0 );
ok( $obj->remove_CytoPositions(), 0 );




my $co1 = Bio::Phenotype::Correlate->new();
my $co2 = Bio::Phenotype::Correlate->new();

ok( $co1->name( "name" ) );

ok( $obj->each_Correlate(), 0 );

ok( $obj->add_Correlates( ( $co1, $co2 ) ) );
ok( $obj->each_Correlate(), 2 );
my @cos = $obj->each_Correlate();
ok( $cos[ 0 ], $co1 );
ok( $cos[ 1 ], $co2 );
ok( $cos[ 0 ]->name, "name" );
ok( $obj->each_Correlate(), 2 );

my @cos2 = $obj->remove_Correlates();
ok( $cos2[ 0 ], $co1 );
ok( $cos2[ 1 ], $co2 );

ok( $obj->each_Correlate(), 0 );
ok( $obj->remove_Correlates(), 0 );




my $m1 = Bio::Phenotype::Measure->new();
my $m2 = Bio::Phenotype::Measure->new();

ok( $m1->description( "desc" ) );

ok( $obj->each_Measure(), 0 );

ok( $obj->add_Measures( ( $m1, $m2 ) ) );
ok( $obj->each_Measure(), 2 );
my @ms = $obj->each_Measure();
ok( $ms[ 0 ], $m1 );
ok( $ms[ 1 ], $m2 );
ok( $ms[ 0 ]->description, "desc" );
ok( $obj->each_Measure(), 2 );

my @ms2 = $obj->remove_Measures();
ok( $ms2[ 0 ], $m1 );
ok( $ms2[ 1 ], $m2 );

ok( $obj->each_Measure(), 0 );
ok( $obj->remove_Measures(), 0 );



ok( $obj->each_keyword(), 0 );

ok( $obj->add_keywords( ( "A", "B" ) ) );
ok( $obj->each_keyword(), 2 );
my @ks = $obj->each_keyword();
ok( $ks[ 0 ], "A" );
ok( $ks[ 1 ], "B" );
ok( $obj->each_keyword(), 2 );

my @ks2 = $obj->remove_keywords();
ok( $ks2[ 0 ], "A" );
ok( $ks2[ 1 ], "B" );

ok( $obj->each_keyword(), 0 );
ok( $obj->remove_keywords(), 0 );



my $l1 = Bio::Annotation::DBLink->new();
my $l2 = Bio::Annotation::DBLink->new();

ok( $l1->comment( "comment" ) );

ok( $obj->each_DBLink(), 0 );

ok( $obj->add_DBLinks( ( $l1, $l2 ) ) );
ok( $obj->each_DBLink(), 2 );
my @ls = $obj->each_DBLink();
ok( $ls[ 0 ], $l1 );
ok( $ls[ 1 ], $l2 );
ok( $ls[ 0 ]->comment(), "comment" );
ok( $obj->each_DBLink(), 2 );

my @ls2 = $obj->remove_DBLinks();
ok( $ls2[ 0 ], $l1 );
ok( $ls2[ 1 ], $l2 );

ok( $obj->each_DBLink(), 0 );
ok( $obj->remove_DBLinks(), 0 );



ok( $obj->each_Genotype(), 0 );

ok( $obj->add_Genotypes( ( "A", "B" ) ) );
ok( $obj->each_Genotype(), 2 );
my @gts = $obj->each_Genotype();
ok( $gts[ 0 ], "A" );
ok( $gts[ 1 ], "B" );
ok( $obj->each_Genotype(), 2 );

my @gts2 = $obj->remove_Genotypes();
ok( $gts2[ 0 ], "A" );
ok( $gts2[ 1 ], "B" );

ok( $obj->each_Genotype(), 0 );
ok( $obj->remove_Genotypes(), 0 );













