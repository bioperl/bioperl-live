# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 12,
               -requires_module => 'Graph::Directed');
	
    use_ok('Bio::Ontology::Relationship');
    use_ok('Bio::Ontology::GOterm');
    use_ok('Bio::Ontology::RelationshipType');
}

my $IS_A = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
isa_ok($IS_A, "Bio::Ontology::RelationshipType");

my $parent = Bio::Ontology::GOterm->new();
isa_ok($parent, "Bio::Ontology::GOterm" );

my $child = Bio::Ontology::GOterm->new();
isa_ok($child, "Bio::Ontology::GOterm");

$parent->name( "parent" );

$child->name( "child" );

my $rel = Bio::Ontology::Relationship->new( -identifier        => "16847",
                                            -parent_term       => $parent,
                                            -child_term        => $child,
                                            -relationship_type => $IS_A );

isa_ok($rel, "Bio::Ontology::Relationship");

is( $rel->identifier(), "16847" );

is( $rel->parent_term()->name(), "parent" );

is( $rel->child_term()->name(), "child" );

is( $rel->relationship_type()->name(), "IS_A" );

ok( $rel->to_string() );
