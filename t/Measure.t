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
    plan tests => 20;
}

use Bio::Phenotype::Measure;
  
my $measure = Bio::Phenotype::Measure->new( -context     => "height",
                                            -description => "desc",
                                            -start       => 10,
                                            -end         => 150,
                                            -unit        => "cm",
                                            -comment     => "comment" );

ok( $measure->isa( "Bio::Phenotype::Measure" ) );

ok( $measure->to_string() );

ok( $measure->context(), "height" );
ok( $measure->description(), "desc" );
ok( $measure->start(), 10 );
ok( $measure->end(), 150 );
ok( $measure->unit(), "cm" );
ok( $measure->comment(), "comment" );

$measure->init();

ok( $measure->context(), "" );
ok( $measure->description(), "" );
ok( $measure->start(), "" );
ok( $measure->end(), "" );
ok( $measure->unit(), "" );
ok( $measure->comment(), "" );

ok( $measure->context( "A" ), "A" );
ok( $measure->description( "B" ), "B" );
ok( $measure->start( "C" ), "C" );
ok( $measure->end( "D" ), "D" );
ok( $measure->unit( "E" ), "E" );
ok( $measure->comment( "F" ), "F" );




