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
    plan tests => 22;
}

use Bio::Ontology::InterProTerm;
use Bio::Ontology::Term;
use Bio::Ontology::Relationship;
use Bio::Ontology::RelationshipType;
use Bio::Ontology::SimpleOntologyEngine;
use Bio::Ontology::Ontology;

my $soe = Bio::Ontology::SimpleOntologyEngine->new();
my $ont = Bio::Ontology::Ontology->new(-name => "Domain",
				       -engine => $soe);

ok( $soe->isa( "Bio::Ontology::OntologyEngineI" ) );
ok ($ont->engine, $soe);

my $ip_term1 = Bio::Ontology::InterProTerm->new( -interpro_id => "IPR000001",
						 -name => "Kringle",
						 -definition => "Kringles are autonomous structural domains ...",
						 -ontology => $ont
					       );

my $ip_term2 = Bio::Ontology::InterProTerm->new( -interpro_id => "IPR000002",
						 -name => "Cdc20/Fizzy",
						 -definition => "The Cdc20/Fizzy region is almost always ...",
						 -ontology => $ont
					       );

my $ip_term3 = Bio::Ontology::InterProTerm->new( -interpro_id => "IPR000003",
						 -name => "Retinoid X receptor",
						 -definition => "Steroid or nuclear hormone receptors ...",
						 -ontology => "Family"
					       );
my $ip_term4 = Bio::Ontology::InterProTerm->new( -interpro_id => "IPR000004",
						 -name => "Test4",
						 -definition => "Test4 definition ...",
						 -ontology => $ont
					       );

# print $ip_term3->name."\n";

my $rel_type = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
my $rel_type1 = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );

my $rel = Bio::Ontology::Relationship->new( -parent_term => $ip_term1,
					    -child_term => $ip_term2,
					    -relationship_type => $rel_type
					  );
my $rel1 = Bio::Ontology::Relationship->new( -parent_term => $ip_term2,
					     -child_term => $ip_term3,
					     -relationship_type => $rel_type
					   );
my $rel2 = Bio::Ontology::Relationship->new( -parent_term => $ip_term1,
					     -child_term => $ip_term4,
					     -relationship_type => $rel_type
					   );
my $rel3 = Bio::Ontology::Relationship->new( -parent_term => $ip_term4,
					     -child_term => $ip_term3,
					     -relationship_type => $rel_type
					   );

$ont->add_term($ip_term1);
$ont->add_term($ip_term2);
$ont->add_term($ip_term3);
$ont->add_term($ip_term4);
$ont->add_relationship($rel);
$ont->add_relationship($rel1);
$ont->add_relationship($rel2);
$ont->add_relationship($rel3);

my @child_terms = $ont->get_child_terms($ip_term1);
ok (scalar(@child_terms), 2);
ok( $child_terms[0], $ip_term2 );
my @child_terms1 = $ont->get_child_terms($ip_term1, $rel_type);
ok (scalar(@child_terms), 2);
ok( $child_terms1[0], $ip_term2 );
ok (scalar($ont->get_child_terms($ip_term1, $rel_type1)), 0);

my @descendant_terms = $ont->get_descendant_terms($ip_term1);
ok( scalar(@descendant_terms), 3);
ok( $descendant_terms[1], $ip_term3 );

my @descendant_terms1 = $ont->get_descendant_terms($ip_term1, $rel_type);
ok( $descendant_terms1[1], $ip_term3 );
ok (scalar(@descendant_terms1), 3);
ok (scalar($ont->get_descendant_terms($ip_term1, $rel_type1)), 0);

my @parent_terms = $ont->get_parent_terms($ip_term2);
ok (scalar(@parent_terms), 1);
ok( $parent_terms[0], $ip_term1 );

my @ancestor_terms = $ont->get_ancestor_terms($ip_term3);
ok( $ancestor_terms[0], $ip_term1 );
ok (scalar(@ancestor_terms), 3);
ok (scalar($ont->get_ancestor_terms($ip_term3, $rel_type)), 3);
ok (scalar($ont->get_ancestor_terms($ip_term3, $rel_type1)), 0);

my @leaf_terms = $ont->get_leaf_terms();
# print scalar(@leaf_terms)."\n";
ok (scalar(@leaf_terms), 1);
ok( $leaf_terms[0], $ip_term3);

my @root_terms = $ont->get_root_terms();
# print scalar(@root_terms)."\n";
ok (scalar(@root_terms), 1);
ok( $root_terms[0], $ip_term1);

#print $ont->engine->to_string();
