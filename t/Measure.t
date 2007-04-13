# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 21;
    use_ok 'Bio::Phenotype::Measure';
}
  
my $measure = Bio::Phenotype::Measure->new( -context     => "height",
                                            -description => "desc",
                                            -start       => 10,
                                            -end         => 150,
                                            -unit        => "cm",
                                            -comment     => "comment" );

isa_ok( $measure, "Bio::Phenotype::Measure" );

ok( $measure->to_string() );

is( $measure->context(), "height" );
is( $measure->description(), "desc" );
is( $measure->start(), 10 );
is( $measure->end(), 150 );
is( $measure->unit(), "cm" );
is( $measure->comment(), "comment" );

$measure->init();

is( $measure->context(), "" );
is( $measure->description(), "" );
is( $measure->start(), "" );
is( $measure->end(), "" );
is( $measure->unit(), "" );
is( $measure->comment(), "" );

is( $measure->context( "A" ), "A" );
is( $measure->description( "B" ), "B" );
is( $measure->start( "C" ), "C" );
is( $measure->end( "D" ), "D" );
is( $measure->unit( "E" ), "E" );
is( $measure->comment( "F" ), "F" );




