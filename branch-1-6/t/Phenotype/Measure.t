# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 21);
	
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
