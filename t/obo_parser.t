# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 45,
			   -requires_module => 'Graph');
	
	use_ok('Bio::OntologyIO');
	use_ok('Bio::Ontology::RelationshipType');
}

my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );

my $parser = Bio::OntologyIO->new(
                      -format    => "obo",
		      -file      => test_input_file('so.obo'));

my $ont = $parser->next_ontology();
ok ($ont);
is ($ont->name(), "sequence");

my @roots = $ont->get_root_terms();
is (scalar(@roots), 1);
is ($roots[0]->name(), "Sequence_Ontology");
is ($roots[0]->identifier(), "SO:0000000");

my @terms = $ont->get_child_terms($roots[0]);
is (scalar(@terms), 5);
my ($term) = grep { $_->name() eq "variation_operation"; } @terms;
ok $term;
($term) = grep { $_->name() eq "sequence_attribute"; } @terms;
ok $term;
($term) = grep { $_->name() eq "consequences_of_mutation"; } @terms;
ok $term;
($term) = grep { $_->name() eq "chromosome_variation"; } @terms;
ok $term;
($term) = grep { $_->name() eq "located_sequence_feature"; } @terms;
ok $term;

@terms = $ont->get_child_terms($terms[0]);
is (scalar(@terms), 5);
($term) = grep { $_->name() eq "translocate"; } @terms;
ok $term;
($term) = grep { $_->name() eq "delete"; } @terms;
ok $term;
($term) = grep { $_->name() eq "insert"; } @terms;
ok $term;
($term) = grep { $_->name() eq "substitute"; } @terms;
ok $term;
($term) = grep { $_->name() eq "invert"; } @terms;
ok $term;

my $featterm = $terms[0];
@terms = $ont->get_child_terms($featterm);
is (scalar(@terms), 2);

# substitution has two parents, see whether this is handled
@terms = $ont->find_terms(-name => "substitution");
$term =  $terms[0];
is ($term->name(), "substitution");

# search using obo terms;
@terms = $ont->find_identically_named_terms($term);
is scalar @terms, 1;
@terms = $ont->find_identical_terms($term);
is scalar @terms, 1;
@terms = $ont->find_similar_terms($term);
is scalar @terms, 7;

@terms = $ont->get_ancestor_terms($term);
is (scalar(@terms), 6);
is (scalar(grep { $_->name() eq "region"; } @terms), 1);
is (scalar(grep { $_->name() eq "sequence_variant"; } @terms), 1);

# processed_transcript has part-of and is-a children

@terms = $ont->find_terms(-name => "processed_transcript");;
$term = $terms[0];

@terms = $ont->get_child_terms($term);
is (scalar(@terms), 5);
@terms = $ont->get_child_terms($term, $PART_OF);
is (scalar(@terms), 2);
@terms = $ont->get_child_terms($term, $IS_A);
is (scalar(@terms), 3);
@terms = $ont->get_child_terms($term, $PART_OF, $IS_A);
is (scalar(@terms), 5);

# TF_binding_site has 2 parents and different relationships in the two
# paths up (although the relationships to its two parents are of the
# same type, namely is-a)
@terms = $ont->find_terms(-name => "TF_binding_site");;
$term = $terms[0];

@terms = $ont->get_parent_terms($term);
is (scalar(@terms), 2);
my ($pterm) = grep { $_->name eq "regulatory_region"; } @terms;
ok $pterm;
@terms = $ont->get_parent_terms($term, $PART_OF);
is (scalar(@terms), 0);
@terms = $ont->get_parent_terms($term, $IS_A);
is (scalar(@terms), 2);
@terms = $ont->get_parent_terms($term, $PART_OF, $IS_A);
is (scalar(@terms), 2);


# pull out all relationships
my @rels = $ont->get_relationships();
my @relset = grep { $_->object_term->name eq "Sequence_Ontology"; } @rels;
is (scalar(@relset), 5);
@relset = grep { $_->subject_term->name eq "Sequence_Ontology"; } @rels;
is (scalar(@relset), 0);

# relationships for a specific term only
($term) = $ont->find_terms(-identifier => "SO:0000082");
ok ($term);
is ($term->identifier, "SO:0000082");
is ($term->name, "processed_transcript_attribute");
@rels = $ont->get_relationships($term);
is (scalar(@rels), 5);
@relset = grep { $_->predicate_term->name eq "IS_A"; } @rels;
is (scalar(@relset), 5);
@relset = grep { $_->object_term->identifier eq "SO:0000082"; } @rels;
is (scalar(@relset), 4);
@relset = grep { $_->subject_term->identifier eq "SO:0000082"; } @rels;
is (scalar(@relset), 1);
