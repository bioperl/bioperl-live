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
    plan tests => 8;
}

use Bio::Ontology::Relationship;
use Bio::Ontology::GOterm;  
use Bio::Ontology::RelationshipType;

my $IS_A = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
ok( $IS_A->isa( "Bio::Ontology::RelationshipType" ) );

my $parent = Bio::Ontology::GOterm->new();
ok( $parent->isa( "Bio::Ontology::GOterm" ) );

my $child = Bio::Ontology::GOterm->new();
ok( $child->isa( "Bio::Ontology::GOterm" ) );

$parent->name( "parent" );

$child->name( "child" );


my $rel = Bio::Ontology::Relationship->new( -identifier        => "16847",
                                            -parent_term       => $parent,
                                            -child_term        => $child,
                                            -relationship_type => $IS_A ); 



ok( $rel->identifier(), "16847" );


ok( $rel->parent_term()->name(), "parent" );

ok( $rel->child_term()->name(), "child" );

ok( $rel->relationship_type()->name(), "IS_A" );

ok( $rel->to_string() );



