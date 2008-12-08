# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 116);
	
    use_ok('Bio::Phenotype::Phenotype');
    use_ok('Bio::Species');
    use_ok('Bio::Annotation::Reference');
    use_ok('Bio::Map::CytoPosition');
    use_ok('Bio::Phenotype::Correlate');
    use_ok('Bio::Phenotype::Measure');
    use_ok('Bio::Annotation::DBLink');
}

my $obj = Bio::Phenotype::Phenotype->new();

isa_ok( $obj, "Bio::Phenotype::PhenotypeI" );
isa_ok( $obj, "Bio::Phenotype::Phenotype" );

ok( $obj->name( "r1" ) );
is( $obj->name(), "r1" );

ok( $obj->description( "This is ..." ) );
is( $obj->description(), "This is ..." );

my $mouse = Bio::Species->new();
$mouse->classification( qw( musculus Mus ) );
ok( $obj->species( $mouse ) );
is( $obj->species()->binomial(), "Mus musculus" );

ok( $obj->comment( "putative" ) );
is( $obj->comment(), "putative" );



is( $obj->each_gene_symbol(), 0 );

ok( $obj->add_gene_symbols( ( "A", "B" ) ) );
is( $obj->each_gene_symbol(), 2 );
my @gs = $obj->each_gene_symbol();
is( $gs[ 0 ], "A" );
is( $gs[ 1 ], "B" );
is( $obj->each_gene_symbol(), 2 );

my @gs2 = $obj->remove_gene_symbols();
is( $gs2[ 0 ], "A" );
is( $gs2[ 1 ], "B" );

is( $obj->each_gene_symbol(), 0 );
is( $obj->remove_gene_symbols(), 0 );



my $v1 = Bio::Variation::VariantI->new();
my $v2 = Bio::Variation::VariantI->new();

$v1->length( "123" );

is( $obj->each_Variant(), 0 );

ok( $obj->add_Variants( ( $v1, $v2 ) ) );
is( $obj->each_Variant(), 2 );
my @vs = $obj->each_Variant();
is( $vs[ 0 ], $v1 );
is( $vs[ 1 ], $v2 );
is( $vs[ 0 ]->length(), "123" );
is( $obj->each_Variant(), 2 );

my @vs2 = $obj->remove_Variants();
is( $vs2[ 0 ], $v1 );
is( $vs2[ 1 ], $v2 );

is( $obj->each_Variant(), 0 );
is( $obj->remove_Variants(), 0 );




my $r1 = Bio::Annotation::Reference->new();
my $r2 = Bio::Annotation::Reference->new();

$r1->title( "title" );

is( $obj->each_Reference(), 0 );

ok( $obj->add_References( ( $r1, $r2 ) ) );
is( $obj->each_Reference(), 2 );
my @rs = $obj->each_Reference();
is( $rs[ 0 ]->display_text, $r1->display_text,'operator overloading in AnnotationI is deprecated');
is( $rs[ 1 ]->display_text, $r2->display_text,'operator overloading in AnnotationI is deprecated');
is( $rs[ 0 ]->title(), "title" );
is( $obj->each_Reference(), 2 );

my @rs2 = $obj->remove_References();
is( $rs2[ 0 ]->display_text, $r1->display_text,'operator overloading in AnnotationI is deprecated');
is( $rs2[ 1 ]->display_text, $r2->display_text,'operator overloading in AnnotationI is deprecated');

is( $obj->each_Reference(), 0 );
is( $obj->remove_References(), 0 );




my $c1 = Bio::Map::CytoPosition->new();
my $c2 = Bio::Map::CytoPosition->new();

$c1->chr( "12" );

is( $obj->each_CytoPosition(), 0 );

ok( $obj->add_CytoPositions( ( $c1, $c2 ) ) );
is( $obj->each_CytoPosition(), 2 );
my @cs = $obj->each_CytoPosition();
is( $cs[ 0 ], $c1 );
is( $cs[ 1 ], $c2 );
is( $cs[ 0 ]->chr(), 12 );
is( $obj->each_CytoPosition(), 2 );

my @cs2 = $obj->remove_CytoPositions();
is( $cs2[ 0 ], $c1 );
is( $cs2[ 1 ], $c2 );

is( $obj->each_CytoPosition(), 0 );
is( $obj->remove_CytoPositions(), 0 );




my $co1 = Bio::Phenotype::Correlate->new();
my $co2 = Bio::Phenotype::Correlate->new();

ok( $co1->name( "name" ) );

is( $obj->each_Correlate(), 0 );

ok( $obj->add_Correlates( ( $co1, $co2 ) ) );
is( $obj->each_Correlate(), 2 );
my @cos = $obj->each_Correlate();
is( $cos[ 0 ], $co1 );
is( $cos[ 1 ], $co2 );
is( $cos[ 0 ]->name, "name" );
is( $obj->each_Correlate(), 2 );

my @cos2 = $obj->remove_Correlates();
is( $cos2[ 0 ], $co1 );
is( $cos2[ 1 ], $co2 );

is( $obj->each_Correlate(), 0 );
is( $obj->remove_Correlates(), 0 );




my $m1 = Bio::Phenotype::Measure->new();
my $m2 = Bio::Phenotype::Measure->new();

ok( $m1->description( "desc" ) );

is( $obj->each_Measure(), 0 );

ok( $obj->add_Measures( ( $m1, $m2 ) ) );
is( $obj->each_Measure(), 2 );
my @ms = $obj->each_Measure();
is( $ms[ 0 ], $m1 );
is( $ms[ 1 ], $m2 );
is( $ms[ 0 ]->description, "desc" );
is( $obj->each_Measure(), 2 );

my @ms2 = $obj->remove_Measures();
is( $ms2[ 0 ], $m1 );
is( $ms2[ 1 ], $m2 );

is( $obj->each_Measure(), 0 );
is( $obj->remove_Measures(), 0 );



is( $obj->each_keyword(), 0 );

ok( $obj->add_keywords( ( "A", "B" ) ) );
is( $obj->each_keyword(), 2 );
my @ks = $obj->each_keyword();
is( $ks[ 0 ], "A" );
is( $ks[ 1 ], "B" );
is( $obj->each_keyword(), 2 );

my @ks2 = $obj->remove_keywords();
is( $ks2[ 0 ], "A" );
is( $ks2[ 1 ], "B" );

is( $obj->each_keyword(), 0 );
is( $obj->remove_keywords(), 0 );



my $l1 = Bio::Annotation::DBLink->new();
my $l2 = Bio::Annotation::DBLink->new();

ok( $l1->comment( "comment" ) );

is( $obj->each_DBLink(), 0 );

ok( $obj->add_DBLinks( ( $l1, $l2 ) ) );
is( $obj->each_DBLink(), 2 );
my @ls = $obj->each_DBLink();
is( $ls[ 0 ]->display_text, $l1->display_text,'operator overloading in AnnotationI is deprecated');
is( $ls[ 1 ]->display_text, $l2->display_text,'operator overloading in AnnotationI is deprecated');
is( $ls[ 0 ]->comment(), "comment" );
is( $obj->each_DBLink(), 2 );

my @ls2 = $obj->remove_DBLinks();
is( $ls2[ 0 ]->display_text, $l1->display_text,'operator overloading in AnnotationI is deprecated');
is( $ls2[ 1 ]->display_text, $l2->display_text,'operator overloading in AnnotationI is deprecated');

is( $obj->each_DBLink(), 0 );
is( $obj->remove_DBLinks(), 0 );



is( $obj->each_Genotype(), 0 );

ok( $obj->add_Genotypes( ( "A", "B" ) ) );
is( $obj->each_Genotype(), 2 );
my @gts = $obj->each_Genotype();
is( $gts[ 0 ], "A" );
is( $gts[ 1 ], "B" );
is( $obj->each_Genotype(), 2 );

my @gts2 = $obj->remove_Genotypes();
is( $gts2[ 0 ], "A" );
is( $gts2[ 1 ], "B" );

is( $obj->each_Genotype(), 0 );
is( $obj->remove_Genotypes(), 0 );
