# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($HAVEGRAPHDIRECTED $NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    eval {
	require Graph::Directed; 
	$HAVEGRAPHDIRECTED=1;
    };
    if ($@) {
	$HAVEGRAPHDIRECTED = 0;
    }

    plan tests => ($NUMTESTS = 9);
}
END {
    foreach ( $Test::ntest..$NUMTESTS) {
	skip('Cannot run tests as Graph::Directed is not installed',1);
   }
}
exit(0) unless $HAVEGRAPHDIRECTED;

require Bio::Ontology::Relationship;
require Bio::Ontology::GOterm;  
require Bio::Ontology::RelationshipType;

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

ok( $rel->isa( "Bio::Ontology::Relationship" ) );

ok( $rel->identifier(), "16847" );

ok( $rel->parent_term()->name(), "parent" );

ok( $rel->child_term()->name(), "child" );

ok( $rel->relationship_type()->name(), "IS_A" );

ok( $rel->to_string() );



