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

use Bio::Ontology::Term;
use Bio::Ontology::Relationship;
use Bio::Ontology::RelationshipType;
use Bio::Ontology::SimpleOntologyEngine;
use Bio::Ontology::Ontology;

my $ont = Bio::Ontology::Ontology->new(-name => "My Ontology");

my $eng = Bio::Ontology::SimpleOntologyEngine->new();
$ont->engine($eng);
ok( $eng->isa( "Bio::Ontology::OntologyEngineI" ) );
ok ($ont->engine, $eng);

my @terms = (
	     [-identifier => "IPR000001",
	      -name => "Kringle",
	      -definition => "Kringles are autonomous structural domains ...",
	      -ontology => $ont
	      ],
	     [-identifier => "IPR000002",
	      -name => "Cdc20/Fizzy",
	      -definition => "The Cdc20/Fizzy region is almost always ...",
	      -ontology => $ont
	      ],
	     [-identifier => "IPR000003",
	      -name => "Retinoid X receptor",
	      -definition => "Steroid or nuclear hormone receptors ...",
	      -ontology => $ont
	      ],
	     [-identifier => "IPR000004",
	      -name => "Test4",
	      -definition => "Test4 definition ...",
	      -ontology => $ont
	      ],
	     );

for(my $i = 0; $i < @terms; $i++) {
    $terms[$i] = Bio::Ontology::Term->new(@{$terms[$i]});
    $ont->add_term($terms[$i]);
}

my $rel_type = Bio::Ontology::RelationshipType->get_instance("IS_A", $ont);
my $rel_type1 = Bio::Ontology::RelationshipType->get_instance("PART_OF", $ont);

my @rels = (
	    [-parent_term => $terms[0],
	     -child_term => $terms[1],
	     -relationship_type => $rel_type,
	     -ontology => $ont,
	     ],
	    [-parent_term => $terms[1],
	     -child_term => $terms[2],
	     -relationship_type => $rel_type,
	     -ontology => $ont,
	     ],
	    [-parent_term => $terms[0],
	     -child_term => $terms[3],
	     -relationship_type => $rel_type,
	     -ontology => $ont,
	     ],
	    [-parent_term => $terms[3],
	     -child_term => $terms[2],
	     -relationship_type => $rel_type,
	     -ontology => $ont,
	     ],
	    );

for(my $i = 0; $i < @rels; $i++) {
    $rels[$i] = Bio::Ontology::Relationship->new(@{$rels[$i]});
    $ont->add_relationship($rels[$i]);
}

my @child_terms = $ont->get_child_terms($terms[0]);
ok (scalar(@child_terms), 2);
ok( $child_terms[0], $terms[1] );
my @child_terms1 = $ont->get_child_terms($terms[0], $rel_type);
ok (scalar(@child_terms), 2);
ok( $child_terms1[0], $terms[1] );
ok (scalar($ont->get_child_terms($terms[0], $rel_type1)), 0);

my @descendant_terms = $ont->get_descendant_terms($terms[0]);
ok( scalar(@descendant_terms), 3);
ok( $descendant_terms[1], $terms[2] );

my @descendant_terms1 = $ont->get_descendant_terms($terms[0], $rel_type);
ok( $descendant_terms1[1], $terms[2] );
ok (scalar(@descendant_terms1), 3);
ok (scalar($ont->get_descendant_terms($terms[0], $rel_type1)), 0);

my @parent_terms = $ont->get_parent_terms($terms[1]);
ok (scalar(@parent_terms), 1);
ok( $parent_terms[0], $terms[0] );

my @ancestor_terms = $ont->get_ancestor_terms($terms[2]);
ok( $ancestor_terms[0], $terms[0] );
ok (scalar(@ancestor_terms), 3);
ok (scalar($ont->get_ancestor_terms($terms[2], $rel_type)), 3);
ok (scalar($ont->get_ancestor_terms($terms[2], $rel_type1)), 0);

my @leaf_terms = $ont->get_leaf_terms();
# print scalar(@leaf_terms)."\n";
ok (scalar(@leaf_terms), 1);
ok( $leaf_terms[0], $terms[2]);

my @root_terms = $ont->get_root_terms();
# print scalar(@root_terms)."\n";
ok (scalar(@root_terms), 1);
ok( $root_terms[0], $terms[0]);

#print $ont->engine->to_string();
