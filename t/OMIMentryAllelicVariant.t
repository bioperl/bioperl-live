# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## # $Id$

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
    plan tests => 23;
}

use Bio::Phenotype::OMIM::OMIMentryAllelicVariant;

my $av = Bio::Phenotype::OMIM::OMIMentryAllelicVariant->new( -number      => ".0001",
                                                             -title       => "ALCOHOL INTOLERANCE",
                                                             -symbol      => "ALDH2*2",
                                                             -description => "The ALDH2*2-encoded",
                                                             -aa_ori      => "GLU",
                                                             -aa_mut      => "LYS",
                                                             -position    => 487 );  

ok( $av->isa( "Bio::Phenotype::OMIM::OMIMentryAllelicVariant" ) );

ok( $av->to_string() );

ok( $av->number(), ".0001" );
ok( $av->title(), "ALCOHOL INTOLERANCE" );
ok( $av->symbol(), "ALDH2*2" );
ok( $av->description(), "The ALDH2*2-encoded" );
ok( $av->aa_ori(), "GLU" );
ok( $av->aa_mut(), "LYS" );
ok( $av->position(), 487 );

$av->init();

ok( $av->number(), "" );
ok( $av->title(), "" );
ok( $av->symbol(), "" );
ok( $av->description(), "" );
ok( $av->aa_ori(), "" );
ok( $av->aa_mut(), "" );
ok( $av->position(), "" );

ok( $av->number( "A" ), "A" );
ok( $av->title( "B" ), "B" );
ok( $av->symbol( "C" ), "C" );
ok( $av->description( "D" ), "D" );
ok( $av->aa_ori( "E" ), "E" );
ok( $av->aa_mut( "F" ), "F" );
ok( $av->position( "G" ), "G" );




