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

    eval { require 'Graph.pm' };
    if( $@ ) {
	    print STDERR "\nShould have Graph.pm installed...\n\n";
	    plan tests => 1;
	    ok( 1 );
	    exit( 0 );
    }

    plan tests => 81;
}


use Bio::Ontology::simpleGOparser;


my $parser = Bio::Ontology::simpleGOparser->new( -go_defs_file_name    => "./t/data/GO.defs.test",
                                                 -components_file_name => "./t/data/component.ontology.test" );


my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );



my $engine = $parser->parse();



my $term = $engine->get_term( "GO:0018897" );

ok( ( $term->each_dblink() )[ 0 ], "MetaCyc:PWY-681" );
ok( ( $term->each_dblink() )[ 1 ], "UM-BBD_pathwayID:dbt" );
ok( ( $term->each_synonym() )[ 0 ], "murein sacculus" );
ok( ( $term->each_synonym() )[ 1 ], "peptidoglycan" );
ok( $term->category()->name(), "components ontology" );
ok( $term->name(), "dibenzothiophene desulfurization" );





$term = $engine->get_term( "GO:0004796" );

ok( ( $term->each_dblink() )[ 0 ], "EC:5.3.99.5" );
ok( ( $term->each_synonym() )[ 0 ], "cytochrome P450 CYP5" );
ok( ( $term->each_secondary_GO_id() )[ 0 ], "GO:0008400" );
ok( $term->category()->name(), "components ontology" );
ok( $term->name(), "thromboxane-A synthase" );





my @parents = sort( $engine->get_parent_terms( $term ) );

ok( @parents == 2 );

ok( $parents[ 0 ]->GO_id(), "GO:0018895" );
ok( $parents[ 1 ]->GO_id(), "GO:0015034" );


@parents = sort( $engine->get_parent_terms( $term, $PART_OF, $IS_A) );

ok( @parents == 2 );

ok( $parents[ 0 ]->GO_id(), "GO:0018895" );
ok( $parents[ 1 ]->GO_id(), "GO:0015034" );








@parents = sort( $engine->get_parent_terms( "GO:0004796", $IS_A ) );

ok( @parents == 2 );

ok( $parents[ 0 ]->GO_id(), "GO:0018895" );
ok( $parents[ 1 ]->GO_id(), "GO:0015034" );


@parents = sort( $engine->get_parent_terms( "GO:0004796", $PART_OF ) );

ok( @parents == 0 );





my @anc = sort( $engine->get_ancestor_terms( $term ) );

ok( @anc == 3 );

ok( $anc[ 0 ]->GO_id(), "GO:0003673" );
ok( $anc[ 1 ]->GO_id(), "GO:0018895" );
ok( $anc[ 2 ]->GO_id(), "GO:0015034" );


@anc = sort( $engine->get_ancestor_terms( "GO:0004796", $IS_A ) );

ok( @anc == 3 );

ok( $anc[ 0 ]->GO_id(), "GO:0003673" );
ok( $anc[ 1 ]->GO_id(), "GO:0018895" );
ok( $anc[ 2 ]->GO_id(), "GO:0015034" );



@anc = sort( $engine->get_ancestor_terms( "GO:0000666" ) );

ok( @anc == 12 );



@anc = sort( $engine->get_ancestor_terms( "GO:0000666", $IS_A ) );

ok( @anc == 2 );

ok( $anc[ 0 ]->GO_id(), "GO:0005811" );
ok( $anc[ 1 ]->GO_id(), "GO:0030481" );





@anc = sort( $engine->get_ancestor_terms( "GO:0000666", $PART_OF ) );

ok( @anc == 6 );
#foreach my $x ( @anc ) {
#    print $x->GO_id(), "\n";
#}

ok( $anc[ 0 ]->GO_id(), "GO:0005623" );
ok( $anc[ 1 ]->GO_id(), "GO:0005933" );
ok( $anc[ 2 ]->GO_id(), "GO:0005935" );
ok( $anc[ 3 ]->GO_id(), "GO:0005938" );
ok( $anc[ 4 ]->GO_id(), "GO:0005937" );
ok( $anc[ 5 ]->GO_id(), "GO:0005625" );





my @childs = sort( $engine->get_child_terms( "GO:0005625", $PART_OF ) );

ok( @childs == 2 );

ok( $childs[ 0 ]->GO_id(), "GO:0000666" );
ok( $childs[ 0 ]->name(), "polarisomeX" );
ok( $childs[ 1 ]->GO_id(), "GO:0000667" );
ok( $childs[ 1 ]->name(), "polarisomeY" );
ok( $childs[ 1 ]->category()->name(), "components ontology" );



ok( $engine->get_term( "GO:0005625" )->name(), "soluble fraction" ); 



@childs = sort( $engine->get_descendant_terms( "GO:0005624", $IS_A ) );



ok( @childs == 6 );

ok( $childs[ 0 ]->GO_id(), "GO:0005792" );
ok( $childs[ 0 ]->name(), "microsome" );
ok( $childs[ 1 ]->GO_id(), "GO:0000299" );
ok( $childs[ 1 ]->name(), "integral membrane protein of membrane fraction" );
ok( $childs[ 2 ]->GO_id(), "GO:0019718" );
ok( $childs[ 2 ]->name(), "rough microsome" );
ok( $childs[ 3 ]->GO_id(), "GO:0019719" );
ok( $childs[ 3 ]->name(), "smooth microsome" );
ok( $childs[ 4 ]->GO_id(), "GO:0000300" );
ok( $childs[ 4 ]->name(), "peripheral membrane protein of membrane fraction" );
ok( $childs[ 5 ]->GO_id(), "GO:0019717" );
ok( $childs[ 5 ]->name(), "synaptosome" );





@childs = sort( $engine->get_descendant_terms( "GO:0005625", $IS_A ) );

ok( @childs == 0 );


@childs = sort( $engine->get_descendant_terms( "GO:0005625", $PART_OF ) );

ok( @childs == 2 );




my @rels = $engine->get_relationships( "GO:0005625" );

ok( @rels == 3 );

ok( $rels[ 0 ]->parent_term()->GO_id(), "GO:0005625" );
ok( $rels[ 0 ]->child_term()->GO_id(), "GO:0000666" );
ok( $rels[ 0 ]->relationship_type()->equals( $PART_OF ) );

ok( $rels[ 1 ]->parent_term()->GO_id(), "GO:0005625" );
ok( $rels[ 1 ]->child_term()->GO_id(), "GO:0000667" );
ok( $rels[ 1 ]->relationship_type()->equals( $PART_OF ) );

ok( $rels[ 2 ]->parent_term()->GO_id(), "GO:0000267" );
ok( $rels[ 2 ]->child_term()->GO_id(), "GO:0005625" );
ok( $rels[ 2 ]->relationship_type()->equals( $IS_A ) );



ok( $engine->graph() );

ok( $engine->add_term( Bio::Ontology::GOterm->new() ) );

ok( $engine->has_term( "GO:0000300" ) );


ok( scalar $engine->each_term(), "44" );

ok( scalar $engine->get_relationship_types(), "2" );

ok( ! $engine->add_relationship( $rels[ 2 ] ) ); # this edge already exists, cannot add

$rels[ 2 ]->child_term()->GO_id( "GO:0005938" );

ok( $engine->add_relationship( $rels[ 2 ] ) ); # now it's changed, can add
 



$parser = Bio::Ontology::simpleGOparser->new( -go_defs_file_name    => "/home/czmasek/GO/GO.defs.test2",
                                              -components_file_name => "/home/czmasek/GO/component.ontology.test2" );



$engine = $parser->parse();

my @roots = sort( $engine->get_root_terms() );

ok( @roots == 1 );

my @leafs = sort( $engine->get_leaf_terms() );

ok( @leafs == 4 );



