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

    plan tests => 50;
}


use Bio::OntologyIO;
use Bio::Ontology::RelationshipType;
use Bio::Root::IO;

my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );

my $io = Bio::Root::IO->new(); # less typing from now on
my $parser = Bio::OntologyIO->new(
                      -format    => "soflat",
		      -file      => $io->catfile("t", "data",
						 "sofa.ontology"),
		      # test overwriting the default
		      -ontology_name => "Sequence Feature Ontology");

my $ont = $parser->next_ontology();
ok ($ont);
ok ($ont->name, "Sequence Feature Ontology");

my @roots = $ont->get_root_terms();
ok (scalar(@roots), 1);
ok ($roots[0]->name(), "Sequence_Feature_Ontology");
ok ($roots[0]->identifier(), "SO:0000000");

my @terms = $ont->get_child_terms($roots[0]);
ok (scalar(@terms), 1);
ok ($terms[0]->name(), "sofa");
@terms = $ont->get_child_terms($terms[0]);
ok (scalar(@terms), 1);
ok ($terms[0]->name(), "feature");
my $featterm = $terms[0];
@terms = $ont->get_child_terms($featterm);
ok (scalar(@terms), 10);

# oligonucleotide has two parents, see whether this is handled
@terms = $ont->get_descendant_terms($featterm);
my ($term) = grep { $_->name() eq "oligonucleotide"; } @terms;
ok $term;
skip(! $term, $term->identifier(), "SO:0000696");

@terms = $ont->get_ancestor_terms($term);
ok (scalar(@terms), 7);
ok (scalar(grep { $_->name() eq "remark"; } @terms), 1);
ok (scalar(grep { $_->name() eq "reagent"; } @terms), 1);

# processed_transcript has part-of and is-a children
@terms = $ont->get_descendant_terms($featterm);
($term) = grep { $_->name() eq "processed_transcript"; } @terms;
ok $term;
skip(! $term, $term->identifier(), "SO:0000233");

@terms = $ont->get_child_terms($term);
ok (scalar(@terms), 4);
@terms = $ont->get_child_terms($term, $PART_OF);
ok (scalar(@terms), 2);
@terms = $ont->get_child_terms($term, $IS_A);
ok (scalar(@terms), 2);
@terms = $ont->get_child_terms($term, $PART_OF, $IS_A);
ok (scalar(@terms), 4);

# now all descendants:
@terms = $ont->get_descendant_terms($term);
ok (scalar(@terms), 13);
@terms = $ont->get_descendant_terms($term, $PART_OF);
ok (scalar(@terms), 2);
@terms = $ont->get_descendant_terms($term, $IS_A);
ok (scalar(@terms), 5);
@terms = $ont->get_descendant_terms($term, $PART_OF, $IS_A);
ok (scalar(@terms), 13);

# TF_binding_site has 2 parents and different relationships in the two
# paths up (although the relationships to its two parents are of the
# same type, namely is-a)
@terms = $ont->get_descendant_terms($featterm);
($term) = grep { $_->name() eq "TF_binding_site"; } @terms;
ok $term;
skip(! $term, $term->identifier(), "SO:0000235");

@terms = $ont->get_parent_terms($term);
ok (scalar(@terms), 2);
my ($pterm) = grep { $_->name eq "regulatory_region"; } @terms;
ok $pterm;
@terms = $ont->get_parent_terms($term, $PART_OF);
ok (scalar(@terms), 0);
@terms = $ont->get_parent_terms($term, $IS_A);
ok (scalar(@terms), 2);
@terms = $ont->get_parent_terms($term, $PART_OF, $IS_A);
ok (scalar(@terms), 2);

# now all ancestors:
@terms = $ont->get_ancestor_terms($term);
ok (scalar(@terms), 6);
@terms = $ont->get_ancestor_terms($term, $PART_OF);
ok (scalar(@terms), 0);
@terms = $ont->get_ancestor_terms($pterm, $PART_OF);
ok (scalar(@terms), 1);
@terms = $ont->get_ancestor_terms($term, $IS_A);
ok (scalar(@terms), 5);
@terms = $ont->get_ancestor_terms($pterm, $IS_A);
ok (scalar(@terms), 0);
@terms = $ont->get_ancestor_terms($term, $PART_OF, $IS_A);
ok (scalar(@terms), 6);

# pull out all relationships
my @rels = $ont->get_relationships();
my @relset = grep { $_->object_term->name eq "sofa"; } @rels;
ok (scalar(@relset), 1);
@relset = grep { $_->subject_term->name eq "sofa"; } @rels;
ok (scalar(@relset), 1);
@relset = grep { $_->object_term->name eq "feature"; } @rels;
ok (scalar(@relset), 10);
@relset = grep { $_->subject_term->name eq "feature"; } @rels;
ok (scalar(@relset), 1);
@relset = grep { $_->object_term->identifier eq "SO:0000233"; } @rels;
ok (scalar(@relset), 4);
@relset = grep { $_->predicate_term->name eq "IS_A" } @relset;
ok (scalar(@relset), 2);

# relationships for a specific term only
($term) = $ont->find_terms(-identifier => "SO:0000233");
ok ($term);
ok ($term->identifier, "SO:0000233");
ok ($term->name, "processed_transcript");
@rels = $ont->get_relationships($term);
ok (scalar(@rels), 5);
@relset = grep { $_->predicate_term->name eq "IS_A"; } @rels;
ok (scalar(@relset), 3);
@relset = grep { $_->object_term->identifier eq "SO:0000233"; } @rels;
ok (scalar(@relset), 4);
