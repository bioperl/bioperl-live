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
	    print STDERR "\nGraph.pm doesn't seem to be installed on this system -- the GO Parser needs it...\n\n";
	    plan tests => 1;
	    ok( 1 );
	    exit( 0 );
    }

    plan tests => 88;
}


use Bio::OntologyIO::simpleGOparser;
use Bio::Root::IO;

my $io = new Bio::Root::IO; # less typing from now on ...
my $parser = Bio::OntologyIO::simpleGOparser->new(
		      -go_defs_file_name => $io->catfile( "t","data",
                                                  "GO.defs.test" ),
		      -components_file_name => $io->catfile( "t","data",
                                                     "component.ontology.test" ) );


my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );


my @onts = ();
while(my $ont = $parser->next_ontology()) {
    push(@onts, $ont);
}
ok (scalar(@onts), 1);
my $ont = $onts[0];
ok ($ont->isa("Bio::Ontology::OntologyI"));
ok ($ont->name(), "GO Components Ontology");

my $engine = $ont->engine();
ok ($engine->isa("Bio::Ontology::OntologyEngineI"));

my $term = $engine->get_term( "GO:0018897" );

my @dblinks = sort ( $term->get_dblinks() );
my @synos = sort ( $term->get_synonyms() );

ok( $dblinks[ 0 ], "MetaCyc:PWY-681" );
ok( $dblinks[ 1 ], "UM-BBD_pathwayID:dbt" );
ok( $synos[ 0 ], "murein sacculus" );
ok( $synos[ 1 ], "peptidoglycan" );
ok( $term->ontology()->name(), "GO Components Ontology" );
ok( $term->name(), "dibenzothiophene desulfurization" );


@dblinks = ();
@synos = ();
$term = "";

$term = $engine->get_term( "GO:0004796" );

@dblinks = sort ( $term->get_dblinks() );
@synos = sort ( $term->get_synonyms() );
my @sec = sort ( $term->get_secondary_GO_ids() ); 

ok( $dblinks[ 0 ], "EC:5.3.99.5" );
ok( $synos[ 0 ], "cytochrome P450 CYP5" );
ok( $sec[ 0 ], "GO:0008400" );
ok( $term->ontology()->name(), "GO Components Ontology" );
ok( $term->name(), "thromboxane-A synthase" );



my @parents = sort goid ( $ont->get_parent_terms( $term ) );

ok( @parents == 2 );

ok( $parents[ 0 ]->GO_id(), "GO:0015034" );
ok( $parents[ 1 ]->GO_id(), "GO:0018895" );



@parents = ();

@parents = sort goid ( $ont->get_parent_terms( $term, $PART_OF, $IS_A) );

ok( @parents == 2 );

ok( $parents[ 0 ]->GO_id(), "GO:0015034" );
ok( $parents[ 1 ]->GO_id(), "GO:0018895" );




@parents = ();

@parents = sort goid ( $ont->get_parent_terms( "GO:0004796", $IS_A ) );

ok( @parents == 2 );

ok( $parents[ 0 ]->GO_id(), "GO:0015034" );
ok( $parents[ 1 ]->GO_id(), "GO:0018895" );



@parents = ();

@parents = sort goid ( $ont->get_parent_terms( "GO:0004796", $PART_OF ) );

ok( @parents == 0 );





my @anc = sort goid ( $ont->get_ancestor_terms( $term ) );

ok( @anc == 3 );

ok( $anc[ 0 ]->GO_id(), "GO:0003673" );
ok( $anc[ 1 ]->GO_id(), "GO:0015034" );
ok( $anc[ 2 ]->GO_id(), "GO:0018895" );


@anc = ();

@anc = sort goid ( $ont->get_ancestor_terms( "GO:0004796", $IS_A ) );

ok( @anc == 3 );

ok( $anc[ 0 ]->GO_id(), "GO:0003673" );
ok( $anc[ 1 ]->GO_id(), "GO:0015034" );
ok( $anc[ 2 ]->GO_id(), "GO:0018895" );


@anc = sort goid ( $ont->get_ancestor_terms( "GO:0000666" ) );

ok( @anc == 12 );



@anc = sort goid ( $ont->get_ancestor_terms( "GO:0000666", $IS_A ) );

ok( @anc == 2 );

ok( $anc[ 0 ]->GO_id(), "GO:0005811" );
ok( $anc[ 1 ]->GO_id(), "GO:0030481" );






@anc = sort goid ( $ont->get_ancestor_terms( "GO:0000666", $PART_OF ) );

ok( @anc == 6 );

ok( $anc[ 0 ]->GO_id(), "GO:0005623" );
ok( $anc[ 1 ]->GO_id(), "GO:0005625" );
ok( $anc[ 2 ]->GO_id(), "GO:0005933" );
ok( $anc[ 3 ]->GO_id(), "GO:0005935" );
ok( $anc[ 4 ]->GO_id(), "GO:0005937" );
ok( $anc[ 5 ]->GO_id(), "GO:0005938" );




my @childs = sort goid ( $ont->get_child_terms( "GO:0005625", $PART_OF ) );

ok( @childs == 2 );

ok( $childs[ 0 ]->GO_id(), "GO:0000666" );
ok( $childs[ 0 ]->name(), "polarisomeX" );
ok( $childs[ 1 ]->GO_id(), "GO:0000667" );
ok( $childs[ 1 ]->name(), "polarisomeY" );
ok( $childs[ 1 ]->ontology()->name(), "GO Components Ontology" );



ok( $engine->get_term( "GO:0005625" )->name(), "soluble fraction" ); 



@childs = sort goid ( $ont->get_descendant_terms( "GO:0005624", $IS_A ) );



ok( @childs == 6 );

ok( $childs[ 0 ]->GO_id(), "GO:0000299" );
ok( $childs[ 0 ]->name(), "integral membrane protein of membrane fraction" );
ok( $childs[ 1 ]->GO_id(), "GO:0000300" );
ok( $childs[ 1 ]->name(), "peripheral membrane protein of membrane fraction" );
ok( $childs[ 2 ]->GO_id(), "GO:0005792" );
ok( $childs[ 2 ]->name(), "microsome" );
ok( $childs[ 3 ]->GO_id(), "GO:0019717" );
ok( $childs[ 3 ]->name(), "synaptosome" );
ok( $childs[ 4 ]->GO_id(), "GO:0019718" );
ok( $childs[ 4 ]->name(), "rough microsome" );
ok( $childs[ 5 ]->GO_id(), "GO:0019719" );
ok( $childs[ 5 ]->name(), "smooth microsome" );





@childs = sort goid ( $ont->get_descendant_terms( "GO:0005625", $IS_A ) );

ok( @childs == 0 );


@childs = sort goid ( $ont->get_descendant_terms( "GO:0005625", $PART_OF ) );

ok( @childs == 2 );




my @rels = sort child_goid ( $ont->get_relationships( "GO:0005625" ) );

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

ok( $ont->add_term( Bio::Ontology::GOterm->new() ) );

ok( $engine->has_term( "GO:0000300" ) );


ok( scalar $engine->each_term(), "44" );

ok( scalar $ont->get_relationship_types(), "2" );

ok( ! $ont->add_relationship( $rels[ 2 ] ) ); # this edge already exists, cannot add

$rels[ 2 ]->child_term()->GO_id( "GO:0005938" );
ok( $ont->add_relationship( $rels[ 2 ] ) ); # now it's changed, can add
 


my @roots = $ont->get_root_terms();
ok( @roots == 10 );

my @leafs = $ont->get_leaf_terms();
ok( @leafs == 19 );



$parser = Bio::OntologyIO::simpleGOparser->new(
		      -go_defs_file_name => $io->catfile( "t", "data",
                                                  "GO.defs.test2" ),
		      -components_file_name => $io->catfile( "t", "data",
                                                     "component.ontology.test2" ) );

$ont = $parser->next_ontology();
ok ($ont);

my @roots2 = $ont->get_root_terms();
ok( @roots2 == 1 );

my @leafs2 = $ont->get_leaf_terms();
ok( @leafs2 == 4 );



sub goid { num ( $a->GO_id() ) <=> num ( $b->GO_id() ) }

sub child_goid { num ( $a->child_term()->GO_id() ) <=> num ( $b->child_term()->GO_id() ) }

sub num {
    my $x = shift( @_ );
    $x =~ s/\D+//g;
    return $x;
}

