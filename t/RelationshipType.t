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

use Bio::Ontology::RelationshipType;
  
my $IS_A     = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
my $PART_OF  = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );
my $CONTAINS = Bio::Ontology::RelationshipType->get_instance( "CONTAINS" );
my $FOUND_IN = Bio::Ontology::RelationshipType->get_instance( "FOUND_IN" );
my $IS_A2    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );

ok( $IS_A->isa( "Bio::Ontology::RelationshipType" ) );


ok( ! $IS_A->equals( $PART_OF ) );

ok( $IS_A->equals( $IS_A2 ) );

ok( $PART_OF->equals( $PART_OF ) );


ok( $IS_A->identifier(), "IS_A" );

ok( $IS_A->name(), "IS_A" );

ok( $IS_A->definition(), "IS_A relationship type" );

ok( $IS_A->category()->name(), "relationship type" );


ok( $PART_OF->identifier(), "PART_OF" );

ok( $PART_OF->name(), "PART_OF" );

ok( $PART_OF->definition(), "PART_OF relationship type" );

ok( $PART_OF->category()->name(), "relationship type" );


ok( $CONTAINS->identifier(), "CONTAINS" );

ok( $CONTAINS->name(), "CONTAINS" );

ok( $CONTAINS->definition(), "CONTAINS relationship type" );

ok( $CONTAINS->category()->name(), "relationship type" );


ok( $FOUND_IN->identifier(), "FOUND_IN" );

ok( $FOUND_IN->name(), "FOUND_IN" );

ok( $FOUND_IN->definition(), "FOUND_IN relationship type" );

ok( $FOUND_IN->category()->name(), "relationship type" );






