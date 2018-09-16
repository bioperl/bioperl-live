# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 102,
			   -requires_module => 'Graph');
	
	use_ok('Bio::OntologyIO');
}

my $parser = Bio::OntologyIO->new(
                      -format    => "go",
		      -defs_file => test_input_file('GO.defs.test'),
                      # test using -file
		      -file      => test_input_file('component.ontology.test'));


my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );


my @onts = ();
while(my $ont = $parser->next_ontology()) {
    push(@onts, $ont);
}
is (scalar(@onts), 1);
my $ont = $onts[0];
isa_ok $ont, "Bio::Ontology::OntologyI";
is ($ont->name(), "Gene Ontology");

my $engine = $ont->engine();
isa_ok $engine, "Bio::Ontology::OntologyEngineI";

my $term = $engine->get_terms( "GO:0018897" );

# note that all dblinks are now Bio::Annotation::DBLink instances and that all
# *dblink* methods related to Bio::Ontology::Term are deprecated; this is due to
# inconsistencies in those Bio::Ontology::Term methods. Use *dbxref* methods
# instead

my @dblinks = sort {$a->display_text cmp $b->display_text} ( $term->get_dbxrefs() );
my @synos = sort ( $term->get_synonyms() );

is( $dblinks[ 0 ]->display_text, "MetaCyc:PWY-681" );
is( $dblinks[ 1 ]->display_text, "UM-BBD_pathwayID:dbt" );
is( $synos[ 0 ], "murein sacculus" );
is( $synos[ 1 ], "peptidoglycan" );
is( $term->ontology()->name(), "Gene Ontology" );
is( $term->name(), "dibenzothiophene desulfurization" );


$term = $engine->get_terms( "GO:0004796" );
@dblinks = sort ( $term->get_dbxrefs() );
@synos = sort ( $term->get_synonyms() );
my @sec = sort ( $term->get_secondary_GO_ids() ); 

is( $dblinks[ 0 ]->display_text, "EC:5.3.99.5" );
is( $synos[ 0 ], "cytochrome P450 CYP5" );
is( $sec[ 0 ], "GO:0008400" );
is( $term->ontology()->name(), "Gene Ontology" );
is( $term->name(), "thromboxane-A synthase" );

my @parents = sort goid ( $ont->get_parent_terms( $term ) );
is( @parents, 2 );
is( $parents[ 0 ]->GO_id(), "GO:0015034" );
is( $parents[ 1 ]->GO_id(), "GO:0018895" );


@parents = sort goid ( $ont->get_parent_terms( $term, $PART_OF, $IS_A) );

is( @parents, 2 );
is( $parents[ 0 ]->GO_id(), "GO:0015034" );
is( $parents[ 1 ]->GO_id(), "GO:0018895" );


@parents = sort goid ( $ont->get_parent_terms( "GO:0004796", $IS_A ) );
is( @parents, 2 );
is( $parents[ 0 ]->GO_id(), "GO:0015034" );
is( $parents[ 1 ]->GO_id(), "GO:0018895" );


@parents = sort goid ( $ont->get_parent_terms( "GO:0004796", $PART_OF ) );
is( scalar(@parents), 0 );
my @anc = sort goid ( $ont->get_ancestor_terms( $term ) );
is( scalar(@anc), 3 );
is( $anc[ 0 ]->GO_id(), "GO:0003673" );
is( $anc[ 1 ]->GO_id(), "GO:0015034" );
is( $anc[ 2 ]->GO_id(), "GO:0018895" );


@anc = sort goid ( $ont->get_ancestor_terms( "GO:0004796", $IS_A ) );
is( scalar(@anc), 3 );
is( $anc[ 0 ]->GO_id(), "GO:0003673" );
is( $anc[ 1 ]->GO_id(), "GO:0015034" );
is( $anc[ 2 ]->GO_id(), "GO:0018895" );


@anc = sort goid ( $ont->get_ancestor_terms( "GO:0000666" ) );
is( @anc, 12 );

@anc = sort goid ( $ont->get_ancestor_terms( "GO:0000666", $IS_A ) );
is( @anc, 2 );
is( $anc[ 0 ]->GO_id(), "GO:0005811" );
is( $anc[ 1 ]->GO_id(), "GO:0030481" );

@anc = sort goid ( $ont->get_ancestor_terms( "GO:0000666", $PART_OF ) );
is( @anc, 6 );
is( $anc[ 0 ]->GO_id(), "GO:0005623" );
is( $anc[ 1 ]->GO_id(), "GO:0005625" );
is( $anc[ 2 ]->GO_id(), "GO:0005933" );
is( $anc[ 3 ]->GO_id(), "GO:0005935" );
is( $anc[ 4 ]->GO_id(), "GO:0005937" );
is( $anc[ 5 ]->GO_id(), "GO:0005938" );


my @childs = sort goid ( $ont->get_child_terms( "GO:0005625", $PART_OF ) );
is( @childs, 2 );
is( $childs[ 0 ]->GO_id(), "GO:0000666" );
is( $childs[ 0 ]->name(), "polarisomeX" );
is( $childs[ 1 ]->GO_id(), "GO:0000667" );
is( $childs[ 1 ]->name(), "polarisomeY" );
is( $childs[ 1 ]->ontology()->name(), "Gene Ontology" );


is( $engine->get_terms( "GO:0005625" )->name(), "soluble fraction" ); 


@childs = sort goid ( $ont->get_descendant_terms( "GO:0005624", $IS_A ) );
is( @childs, 6 );
is( $childs[ 0 ]->GO_id(), "GO:0000299" );
is( $childs[ 0 ]->name(), "integral membrane protein of membrane fraction" );
is( $childs[ 1 ]->GO_id(), "GO:0000300" );
is( $childs[ 1 ]->name(), "peripheral membrane protein of membrane fraction" );
is( $childs[ 2 ]->GO_id(), "GO:0005792" );
is( $childs[ 2 ]->name(), "microsome" );
is( $childs[ 3 ]->GO_id(), "GO:0019717" );
is( $childs[ 3 ]->name(), "synaptosome" );
is( $childs[ 4 ]->GO_id(), "GO:0019718" );
is( $childs[ 4 ]->name(), "rough microsome" );
is( $childs[ 5 ]->GO_id(), "GO:0019719" );
is( $childs[ 5 ]->name(), "smooth microsome" );

@childs = sort goid ( $ont->get_descendant_terms( "GO:0005625", $IS_A ) );
is( @childs, 0 );


@childs = sort goid ( $ont->get_descendant_terms( "GO:0005625", $PART_OF ) );
is( @childs, 2 );

my @rels = sort child_goid ( $ont->get_relationships( "GO:0005625" ) );
is( @rels, 3 );
is( $rels[ 0 ]->object_term()->GO_id(), "GO:0005625" );
is( $rels[ 0 ]->subject_term()->GO_id(), "GO:0000666" );
ok( $rels[ 0 ]->predicate_term()->equals( $PART_OF ) );

is( $rels[ 1 ]->object_term()->GO_id(), "GO:0005625" );
is( $rels[ 1 ]->subject_term()->GO_id(), "GO:0000667" );
ok( $rels[ 1 ]->predicate_term()->equals( $PART_OF ) );

is( $rels[ 2 ]->object_term()->GO_id(), "GO:0000267" );
is( $rels[ 2 ]->subject_term()->GO_id(), "GO:0005625" );
ok( $rels[ 2 ]->predicate_term()->equals( $IS_A ) );

# dbxrefs and synonyms are candidates for being falsely picked up by
# overly promiscuous regular expressions as related terms, so we test for
# that here
my @terms = $engine->get_terms( "EC:5.3.99.5" );
is (scalar(@terms), 0);
@terms = $engine->get_terms("MetaCyc:PWY-681","MetaCyc:PWY");
is (scalar(@terms), 0);
@terms = $engine->get_terms("UM-BBD_pathwayID:dbt","BBD_pathwayID:dbt",
                            "UM-BBD_pathwayID:dbt2","BBD_pathwayID:dbt2");
is (scalar(@terms), 0);


ok( $engine->graph() );

ok( $ont->add_term( Bio::Ontology::GOterm->new(-identifier => "GO:0000000")));

ok( $engine->has_term( "GO:0000300" ) );

is( scalar $ont->get_all_terms(), 44 );
is( scalar $ont->get_relationship_types(), 3 );

ok( ! $ont->add_relationship( $rels[ 2 ] ) ); # this edge already exists, cannot add

$rels[ 2 ]->subject_term()->GO_id( "GO:0005938" );
ok( $ont->add_relationship( $rels[ 2 ] ) ); # now it's changed, can add
 

my @roots = $ont->get_root_terms();
is( scalar(@roots), 10 );

my @leafs = $ont->get_leaf_terms();
is( scalar(@leafs), 19 );



$parser = Bio::OntologyIO->new(
                      -format    => "go",
		      -defs_file => test_input_file('GO.defs.test2'),
		      # test using -files
		      -files     => test_input_file('component.ontology.test2'));

$ont = $parser->next_ontology();
ok ($ont);

@roots = $ont->get_root_terms();
is( scalar(@roots), 1 );

@leafs = $ont->get_leaf_terms();
is( scalar(@leafs), 4 );

$parser = Bio::OntologyIO->new(
                      -format    => "go",
		      -file      => test_input_file('mpath.ontology.test'));

ok($parser);
$ont = $parser->next_ontology;
ok($ont);
$engine = $ont->engine;
ok($engine);
$term = $engine->get_terms( "MPATH:30" );
is($term->identifier,"MPATH:30");
is($term->name,"cystic medial necrosis");
is($term->definition,undef);
is((sort $term->get_synonyms)[0],"erdheim disease");
is($ont->get_parent_terms( $term )->name,"tissue specific degenerative process");
is(scalar($ont->get_root_terms()),2);
@anc = $ont->get_ancestor_terms($term);
is(scalar(@anc),4);

#################################################################
# helper functions
#################################################################

sub goid { snum ( $a->GO_id() ) <=> snum ( $b->GO_id() ) }

sub child_goid { snum ( $a->child_term()->GO_id() ) <=> snum ( $b->child_term()->GO_id() ) }

sub snum {
    my $x = shift( @_ );
    $x =~ s/\D+//g;
    return $x;
}
