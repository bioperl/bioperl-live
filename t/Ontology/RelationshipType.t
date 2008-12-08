# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 23,
               -requires_module => 'Graph::Directed');
	
    use_ok('Bio::Ontology::RelationshipType');
    use_ok('Bio::Ontology::Ontology');
}

my $ont = Bio::Ontology::Ontology->new(-name => "relationship type");
  
my $IS_A     = Bio::Ontology::RelationshipType->get_instance("IS_A", $ont);
my $PART_OF  = Bio::Ontology::RelationshipType->get_instance("PART_OF", $ont);
my $CONTAINS = Bio::Ontology::RelationshipType->get_instance("CONTAINS", $ont);
my $FOUND_IN = Bio::Ontology::RelationshipType->get_instance("FOUND_IN", $ont);
my $IS_A2    = Bio::Ontology::RelationshipType->get_instance("IS_A", $ont);

isa_ok($IS_A, "Bio::Ontology::RelationshipType");
isa_ok($IS_A, "Bio::Ontology::TermI");


ok( ! $IS_A->equals( $PART_OF ) );
ok( $IS_A->equals( $IS_A2 ) );
ok( $PART_OF->equals( $PART_OF ) );


is( $IS_A->identifier(), undef ); # don't make up identifiers
is( $IS_A->name(), "IS_A" );
is( $IS_A->definition(), "IS_A relationship predicate (type)" );
is( $IS_A->ontology()->name(), "relationship type" );


is( $PART_OF->identifier(), undef ); # don't make up identifiers
is( $PART_OF->name(), "PART_OF" );
is( $PART_OF->definition(), "PART_OF relationship predicate (type)" );
is( $PART_OF->ontology()->name(), "relationship type" );


is( $CONTAINS->identifier(), undef ); # don't make up identifiers
is( $CONTAINS->name(), "CONTAINS" );
is( $CONTAINS->definition(), "CONTAINS relationship predicate (type)" );
is( $CONTAINS->ontology()->name(), "relationship type" );


is( $FOUND_IN->identifier(), undef ); # don't make up identifiers
is( $FOUND_IN->name(), "FOUND_IN" );
is( $FOUND_IN->definition(), "FOUND_IN relationship predicate (type)" );
is( $FOUND_IN->ontology()->name(), "relationship type" );
