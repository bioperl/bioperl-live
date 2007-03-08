# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $ Id: Exp $

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 17;
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




