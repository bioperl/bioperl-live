#!/usr/bin/perl -W

use strict;

use Bio::Ontology::simpleGOparser;


my $parser = Bio::Ontology::simpleGOparser->new( -go_defs_file_name    => "GO.defs",
                                                 -components_file_name => "component.ontology",
                                                 -functions_file_name  => "function.ontology",
                                                 -processes_file_name  => "process.ontology" );



my $engine = $parser->parse();



my $IS_A    = $engine->is_a_relationship();
my $PART_OF = $engine->part_of_relationship();


# root terms
# ----------

print scalar( $engine->get_root_terms() ), "\n";

print( ( $engine->get_root_terms() )[ 0 ]->to_string(), "\n" );

print "\n\n";


# parent terms
# ------------

my @parents1 = $engine->get_parent_terms( "GO:0045474", $IS_A );

my @parents2 = $engine->get_parent_terms( $parents1[ 0 ] );

foreach my $p ( @parents2 ) {
    print $p->to_string(), "\n";
}

print "\n\n";



# child terms
# -----------


my @child0 = $engine->get_child_terms( ( $engine->get_root_terms() )[ 0 ] );

foreach my $c ( @child0 ) {
    print $c->to_string(), "\n";
}

print "\n\n";

my @child1 = $engine->get_child_terms( "GO:0008044", $IS_A );

foreach my $cc ( @child1 ) {
    print $cc->to_string(), "\n";
}

print "\n\n";


# ancestor terms
# --------------

my @ancestor1 = $engine->get_ancestor_terms( "GO:0007323", $PART_OF );



# descendant terms
# ----------------

my @descendant1 = $engine->get_descendant_terms( "GO:0007323", $PART_OF );



# relationships
# -------------

my @rel1 = $engine->get_relationships( "GO:0007323", $PART_OF, $IS_A );


# leaf terms
# ----------

my @l1 = $engine->get_leaf_terms();


# get term with id
# ----------------

my $term = $engine->get_term( "GO:0001661" );


# all terms
# ---------

my @terms = $engine->each_term();


exit( 1 );


