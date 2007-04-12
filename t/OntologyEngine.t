# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($HAVEGRAPHDIRECTED $DEBUG $NUMTESTS);
BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	eval {
		require Graph::Directed;
	};
	if ($@) {
		plan skip_all => "Cannot run tests, Graph::Directed not installed\n";
	} else {
		plan tests => ($NUMTESTS = 27);
	}
	use_ok('Bio::Ontology::Term');
	use_ok('Bio::Ontology::Relationship');
	use_ok('Bio::Ontology::RelationshipType');
	use_ok('Bio::Ontology::SimpleOntologyEngine');
	use_ok('Bio::Ontology::Ontology');
}

my $ont = Bio::Ontology::Ontology->new(-name => "My Ontology");

my $eng = Bio::Ontology::SimpleOntologyEngine->new();
$ont->engine($eng);
isa_ok( $eng,"Bio::Ontology::OntologyEngineI" );
is ($ont->engine, $eng);

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
	    [-object_term => $terms[0],
	     -subject_term => $terms[1],
	     -predicate_term => $rel_type,
	     -ontology => $ont,
	     ],
	    [-object_term => $terms[1],
	     -subject_term => $terms[2],
	     -predicate_term => $rel_type,
	     -ontology => $ont,
	     ],
	    [-object_term => $terms[0],
	     -subject_term => $terms[3],
	     -predicate_term => $rel_type,
	     -ontology => $ont,
	     ],
	    [-object_term => $terms[3],
	     -subject_term => $terms[2],
	     -predicate_term => $rel_type,
	     -ontology => $ont,
	     ],
	    );

for(my $i = 0; $i < @rels; $i++) {
    $rels[$i] = Bio::Ontology::Relationship->new(@{$rels[$i]});
    $ont->add_relationship($rels[$i]);
}

my @child_terms = sort { $a->identifier() cmp $b->identifier();
		     } $ont->get_child_terms($terms[0]);
is (scalar(@child_terms), 2);
is( $child_terms[0], $terms[1] );
my @child_terms1 = sort { $a->identifier() cmp $b->identifier();
		      } $ont->get_child_terms($terms[0], $rel_type);
is (scalar(@child_terms), 2);
is( $child_terms1[0], $terms[1] );
is (scalar($ont->get_child_terms($terms[0], $rel_type1)), 0);

my @descendant_terms = sort { $a->identifier() cmp $b->identifier();
			  } $ont->get_descendant_terms($terms[0]);
is( scalar(@descendant_terms), 3);
is( $descendant_terms[1], $terms[2] );

my @descendant_terms1 = sort { $a->identifier() cmp $b->identifier();
			   } $ont->get_descendant_terms($terms[0], $rel_type);
is( $descendant_terms1[1], $terms[2] );
is (scalar(@descendant_terms1), 3);
is (scalar($ont->get_descendant_terms($terms[0], $rel_type1)), 0);

my @parent_terms = sort { $a->identifier() cmp $b->identifier();
		      } $ont->get_parent_terms($terms[1]);
is (scalar(@parent_terms), 1);
is( $parent_terms[0], $terms[0] );

my @ancestor_terms = sort { $a->identifier() cmp $b->identifier();
			} $ont->get_ancestor_terms($terms[2]);
is( $ancestor_terms[0], $terms[0] );
is (scalar(@ancestor_terms), 3);
is (scalar($ont->get_ancestor_terms($terms[2], $rel_type)), 3);
is (scalar($ont->get_ancestor_terms($terms[2], $rel_type1)), 0);

my @leaf_terms = $ont->get_leaf_terms();
# print scalar(@leaf_terms)."\n";
is (scalar(@leaf_terms), 1);
is( $leaf_terms[0], $terms[2]);

my @root_terms = $ont->get_root_terms();
# print scalar(@root_terms)."\n";
is (scalar(@root_terms), 1);
is( $root_terms[0], $terms[0]);

#print $ont->engine->to_string();
