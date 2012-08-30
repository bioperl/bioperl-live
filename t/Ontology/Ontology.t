# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(
        -tests           => 55,
        -requires_module => 'Graph'
    );

    use_ok('Bio::OntologyIO');
    use_ok('Bio::Ontology::RelationshipType');
}

my $IS_A    = Bio::Ontology::RelationshipType->get_instance("IS_A");
my $PART_OF = Bio::Ontology::RelationshipType->get_instance("PART_OF");

my $parser = Bio::OntologyIO->new(
    -format => "soflat",
    -file   => test_input_file('sofa.ontology')
);

my $ont = $parser->next_ontology();
isa_ok( $ont, 'Bio::Ontology::Ontology' );
is( $ont->name, "Sequence Feature Ontology" );

my @roots = $ont->get_root_terms();
is( scalar(@roots),          1 );
is( $roots[0]->name(),       "Sequence_Feature_Ontology" );
is( $roots[0]->identifier(), "SO:0000000" );

my @terms = $ont->get_child_terms( $roots[0] );
is( scalar(@terms),    1 );
is( $terms[0]->name(), "sofa" );
@terms = $ont->get_child_terms( $terms[0] );
is( scalar(@terms),    1 );
is( $terms[0]->name(), "feature" );
my $featterm = $terms[0];
@terms = $ont->get_child_terms($featterm);
is( scalar(@terms), 10 );

# oligonucleotide has two parents, see whether this is handled
@terms = $ont->get_descendant_terms($featterm);
my ($term) = grep { $_->name() eq "oligonucleotide"; } @terms;
ok $term;

#TODO: {
#	local $TODO = '$term->identifier()';
is( $term->identifier(), "SO:0000696" );

#}

@terms = $ont->get_ancestor_terms($term);
is( scalar(@terms), 7 );
is( scalar( grep { $_->name() eq "remark"; } @terms ),  1 );
is( scalar( grep { $_->name() eq "reagent"; } @terms ), 1 );

# processed_transcript has part-of and is-a children
@terms = $ont->get_descendant_terms($featterm);
($term) = grep { $_->name() eq "processed_transcript"; } @terms;
ok $term;

#TODO: {
#	local $TODO = '$term->identifier()';
is( $term->identifier(), "SO:0000233" );

#}

@terms = $ont->get_child_terms($term);
is( scalar(@terms), 4 );
@terms = $ont->get_child_terms( $term, $PART_OF );
is( scalar(@terms), 2 );
@terms = $ont->get_child_terms( $term, $IS_A );
is( scalar(@terms), 2 );
@terms = $ont->get_child_terms( $term, $PART_OF, $IS_A );
is( scalar(@terms), 4 );

# now all descendants:
@terms = $ont->get_descendant_terms($term);
is( scalar(@terms), 13 );
@terms = $ont->get_descendant_terms( $term, $PART_OF );
is( scalar(@terms), 2 );
@terms = $ont->get_descendant_terms( $term, $IS_A );
is( scalar(@terms), 5 );
@terms = $ont->get_descendant_terms( $term, $PART_OF, $IS_A );
is( scalar(@terms), 13 );

# TF_binding_site has 2 parents and different relationships in the two
# paths up (although the relationships to its two parents are of the
# same type, namely is-a)
@terms = $ont->get_descendant_terms($featterm);
($term) = grep { $_->name() eq "TF_binding_site"; } @terms;
ok $term;

#TODO: {
#	local $TODO = '$term->identifier()';
is( $term->identifier(), "SO:0000235" );

#}

@terms = $ont->get_parent_terms($term);
is( scalar(@terms), 2 );
my ($pterm) = grep { $_->name eq "regulatory_region"; } @terms;
ok $pterm;
@terms = $ont->get_parent_terms( $term, $PART_OF );
is( scalar(@terms), 0 );
@terms = $ont->get_parent_terms( $term, $IS_A );
is( scalar(@terms), 2 );
@terms = $ont->get_parent_terms( $term, $PART_OF, $IS_A );
is( scalar(@terms), 2 );

# now all ancestors:
@terms = $ont->get_ancestor_terms($term);
is( scalar(@terms), 6 );
@terms = $ont->get_ancestor_terms( $term, $PART_OF );
is( scalar(@terms), 0 );
@terms = $ont->get_ancestor_terms( $pterm, $PART_OF );
is( scalar(@terms), 1 );
@terms = $ont->get_ancestor_terms( $term, $IS_A );
is( scalar(@terms), 5 );
@terms = $ont->get_ancestor_terms( $pterm, $IS_A );
is( scalar(@terms), 0 );
@terms = $ont->get_ancestor_terms( $term, $PART_OF, $IS_A );
is( scalar(@terms), 6 );

# pull out all relationships
my @rels = $ont->get_relationships();
my @relset = grep { $_->object_term->name eq "sofa"; } @rels;
is( scalar(@relset), 1 );
@relset = grep { $_->subject_term->name eq "sofa"; } @rels;
is( scalar(@relset), 1 );
@relset = grep { $_->object_term->name eq "feature"; } @rels;
is( scalar(@relset), 10 );
@relset = grep { $_->subject_term->name eq "feature"; } @rels;
is( scalar(@relset), 1 );
@relset = grep { $_->object_term->identifier eq "SO:0000233"; } @rels;
is( scalar(@relset), 4 );
@relset = grep { $_->predicate_term->name eq "IS_A" } @relset;
is( scalar(@relset), 2 );

# relationships for a specific term only
($term) = $ont->find_terms( -identifier => "SO:0000233" );
ok($term);
is( $term->identifier, "SO:0000233" );
is( $term->name,       "processed_transcript" );
@rels = $ont->get_relationships($term);
is( scalar(@rels), 5 );
@relset = grep { $_->predicate_term->name eq "IS_A"; } @rels;
is( scalar(@relset), 3 );
@relset = grep { $_->object_term->identifier eq "SO:0000233"; } @rels;
is( scalar(@relset), 4 );

SKIP: {
    test_skip(-tests    => 3, -requires_module => 'XML::Parser::PerlSAX');
    for my $file (qw/interpro.xml interpro_sample.xml interpro_relationship.xml/) {
        my $parser = Bio::OntologyIO->new(
            -format => "interpro",
            -file   => test_input_file($file),
        );
        ok( $parser->next_ontology, "Interpro XML file $file can be parsed" );
    }
}

