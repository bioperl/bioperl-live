# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 17);
    
    use_ok('Bio::Phenotype::Correlate');
    use_ok('Bio::Species');
}

my $mouse = Bio::Species->new();
  
$mouse->classification( qw( musculus Mus ) );

my $co = Bio::Phenotype::Correlate->new( -name        => "4(Tas1r3)",
                                         -description => "mouse correlate of human phenotype MIM 605865",
                                         -species     => $mouse,
                                         -type        => "homolog",
                                         -comment     => "type=homolog is putative" );

isa_ok($co, "Bio::Phenotype::Correlate" );

ok( $co->to_string() );

is( $co->name(), "4(Tas1r3)" );
is( $co->description(), "mouse correlate of human phenotype MIM 605865" );
is( $co->species()->binomial(), "Mus musculus" );
is( $co->type(), "homolog" );
is( $co->comment(), "type=homolog is putative" );

$co->init();

is( $co->name(), "" );
is( $co->description(), "" );
is( $co->type(), "" );
is( $co->comment(), "" );

is( $co->name( "A" ), "A" );
is( $co->description( "B" ), "B" );
is( $co->type( "C" ), "C" );
is( $co->comment( "D" ), "D" );
