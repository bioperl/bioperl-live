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
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 15;
}

use Bio::Phenotype::Correlate;
use Bio::Species;

my $mouse = Bio::Species->new();
  
$mouse->classification( qw( musculus Mus ) );

my $co = Bio::Phenotype::Correlate->new( -name        => "4(Tas1r3)",
                                         -description => "mouse correlate of human phenotype MIM 605865",
                                         -species     => $mouse,
                                         -type        => "homolog",
                                         -comment     => "type=homolog is putative" );

ok( $co->isa( "Bio::Phenotype::Correlate" ) );

ok( $co->to_string() );

ok( $co->name(), "4(Tas1r3)" );
ok( $co->description(), "mouse correlate of human phenotype MIM 605865" );
ok( $co->species()->binomial(), "Mus musculus" );
ok( $co->type(), "homolog" );
ok( $co->comment(), "type=homolog is putative" );

$co->init();

ok( $co->name(), "" );
ok( $co->description(), "" );
ok( $co->type(), "" );
ok( $co->comment(), "" );

ok( $co->name( "A" ), "A" );
ok( $co->description( "B" ), "B" );
ok( $co->type( "C" ), "C" );
ok( $co->comment( "D" ), "D" );




