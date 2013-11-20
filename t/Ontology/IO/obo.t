# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 92,
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

my @terms = sort {$a->name cmp $b->name} $ont->get_child_terms($roots[0]);
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

@terms = sort {$a->name cmp $b->name} $ont->get_child_terms($terms[4]);
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
@terms = sort {$a->name cmp $b->name} $ont->get_child_terms($featterm);
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




#### --- testing obo parsers for regulates relationships
my $parser2 = Bio::OntologyIO->new (
	 -format => 'obo',
	 -file => test_input_file('regulation_test.obo'));

isa_ok($parser2,'Bio::OntologyIO', 'got a ontology IO handler');

my @ontologies;
while (my $ont = $parser2->next_ontology()) {
	 isa_ok($ont,'Bio::Ontology::Ontology','got ontology parser2');
	 isa_ok($ont->engine,'Bio::Ontology::OBOEngine','got OBO engine object');
	 push @ontologies,$ont;
}

my $molont = $ontologies[1];
my $bioont = $ontologies[2];
is($ontologies[0]->name(),'gene_ontology','Gene ontology');
is($bioont->name(),'biological_process','biological process');
is($molont->name(),'molecular_function','molecular function');

my ($broot) = $bioont->get_root_terms();
is($broot->name(),'biological_process','Got root');


my ($mroot) = $molont->get_root_terms();
is($mroot->name(),'molecular_function','Got root');


## -- testing newly introduced relationships
is($ontologies[0]->get_relationship_type('REGULATES')->name,'REGULATES','Got regulates
	 from gene_ontology');
is($ontologies[0]->get_relationship_type('POSITIVELY_REGULATES')->name,'POSITIVELY_REGULATES','Got
	 positively regulates from gene_ontology');
is($bioont->get_relationship_type('REGULATES')->name,'REGULATES','Got
	  regulates from biological_process');
is($bioont->get_relationship_type('POSITIVELY_REGULATES')->name,'POSITIVELY_REGULATES','Got
	 positively regulates from biological_process');


## -- getting relationships for various ontologies
my @onto_pred = sort {$a->name cmp $b->name} $ontologies[0]->get_predicate_terms();
my @bio_pred =  sort {$a->name cmp $b->name} $bioont->get_predicate_terms();
is(scalar @onto_pred,6,'Got predicates for gene_ontology');
is(scalar @bio_pred,2,'Got predicates for biological_process');
is($onto_pred[4]->name(),'REGULATES','Got regulates predicate');
is($bio_pred[0]->name(),'POSITIVELY_REGULATES','Got positively regulates predicate');


my @bio_rel = $bioont->get_relationships();
my @mol_rel = $molont->get_relationships();
is(scalar @bio_rel,11,'Got relationships for biological_process');
is(scalar @mol_rel,2,'Got relationships for molecular_function');
is($mol_rel[0]->predicate_term->name(),'IS_A','Got is a relationship from
	 molecular_function');
## ----


## -- testing the regulates relationships between term1s
my $REG = Bio::Ontology::RelationshipType->get_instance('REGULATES');
my $PREG = Bio::Ontology::RelationshipType->get_instance('POSITIVELY_REGULATES');

my ($term1) = $bioont->find_terms(-identifier => 'GO:0050790');
isa_ok($term1,'Bio::Ontology::Term','Got term object');
is($term1->identifier(),'GO:0050790', 'Got term id');
is($term1->name(),'regulation of catalytic activity', 'Got term name');

my ($parent) = $bioont->get_parent_terms($term1,$REG);
isa_ok($parent,'Bio::Ontology::Term','Got regulated object');
is($parent->identifier(),'GO:0003824','Got regulated term1 id'); 

## -- now testing the other way
my ($child) = $bioont->get_child_terms($parent,$REG); 
isa_ok($child,'Bio::Ontology::Term','Got term1 object');
is($child->identifier(),$term1->identifier(),'Got back the child');

my ($term2) = $bioont->find_terms(-identifier => 'GO:0043085');
isa_ok($term2,'Bio::Ontology::Term','Got term object');
is($term2->identifier(),'GO:0043085', 'Got term id');
is($term2->name(),'positive regulation of catalytic activity', 'Got term name');

my ($parent2) = $bioont->get_parent_terms($term2,$PREG);
isa_ok($parent2,'Bio::Ontology::Term','Got regulated object');
is($parent2->identifier(),'GO:0003824','Got regulated term1 id'); 
is($parent->name(),$parent2->name(),'Got identical regulation');


my ($child2) = $bioont->get_child_terms($parent2,$PREG); 
isa_ok($child2,'Bio::Ontology::Term','Got term1 object');
is($child2->identifier(),$term2->identifier(),'Got back the child');



#### --- testing obo parsers for secondary identifiers
my $parser3 = Bio::OntologyIO->new (
	 -format => 'obo',
	 -file => test_input_file('sp_subset.obo'));

isa_ok($parser3,'Bio::OntologyIO', 'got a ontology IO handler');


my $sp_ont = $parser3->next_ontology();
ok ($sp_ont);
is ($sp_ont->name(), "solanaceae_phenotype");

#the term 'plant size has 4 secondary identifiers
@terms = $sp_ont->find_terms(-name => "plant size");
$term =  $terms[0];
is ($term->name(), "plant size");

my @xrefs = $term->get_secondary_ids();
is(scalar(@xrefs) , 4);

my ($xref) = grep { $_ eq "PATO:0000117"; } @xrefs;
ok $xref;
($xref) = grep { $_ eq "PO:0000003"; } @xrefs;
ok $xref;
($xref) = grep { $_ eq "PO:0007130"; } @xrefs;
ok $xref;
($xref) = grep { $_ eq "TO:0000207"; } @xrefs;
ok $xref;


