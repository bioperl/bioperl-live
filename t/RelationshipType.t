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
    plan tests => 21;
}

use Bio::Ontology::RelationshipType;
use Bio::Ontology::Ontology;

my $ont = Bio::Ontology::Ontology->new(-name => "relationship type");
  
my $IS_A     = Bio::Ontology::RelationshipType->get_instance("IS_A", $ont);
my $PART_OF  = Bio::Ontology::RelationshipType->get_instance("PART_OF", $ont);
my $CONTAINS = Bio::Ontology::RelationshipType->get_instance("CONTAINS", $ont);
my $FOUND_IN = Bio::Ontology::RelationshipType->get_instance("FOUND_IN", $ont);
my $IS_A2    = Bio::Ontology::RelationshipType->get_instance("IS_A", $ont);

ok( $IS_A->isa( "Bio::Ontology::RelationshipType" ) );
ok( $IS_A->isa( "Bio::Ontology::TermI" ) );


ok( ! $IS_A->equals( $PART_OF ) );
ok( $IS_A->equals( $IS_A2 ) );
ok( $PART_OF->equals( $PART_OF ) );


ok( $IS_A->identifier(), undef ); # don't make up identifiers
ok( $IS_A->name(), "IS_A" );
ok( $IS_A->definition(), "IS_A relationship predicate (type)" );
ok( $IS_A->ontology()->name(), "relationship type" );


ok( $PART_OF->identifier(), undef ); # don't make up identifiers
ok( $PART_OF->name(), "PART_OF" );
ok( $PART_OF->definition(), "PART_OF relationship predicate (type)" );
ok( $PART_OF->ontology()->name(), "relationship type" );


ok( $CONTAINS->identifier(), undef ); # don't make up identifiers
ok( $CONTAINS->name(), "CONTAINS" );
ok( $CONTAINS->definition(), "CONTAINS relationship predicate (type)" );
ok( $CONTAINS->ontology()->name(), "relationship type" );


ok( $FOUND_IN->identifier(), undef ); # don't make up identifiers
ok( $FOUND_IN->name(), "FOUND_IN" );
ok( $FOUND_IN->definition(), "FOUND_IN relationship predicate (type)" );
ok( $FOUND_IN->ontology()->name(), "relationship type" );






