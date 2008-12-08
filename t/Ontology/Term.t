# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 54,
			   -requires_module => 'Graph::Directed');
	
	use_ok('Bio::Ontology::Term');
	use_ok('Bio::Ontology::TermFactory');
	use_ok('Bio::Annotation::DBLink');
}

my $obj = Bio::Ontology::Term->new();

isa_ok $obj, "Bio::Ontology::TermI";

is( $obj->identifier( "0003947" ), "0003947" );
is( $obj->identifier(), "0003947" );

is( $obj->name( "N-acetylgalactosaminyltransferase" ), "N-acetylgalactosaminyltransferase" );
is( $obj->name(), "N-acetylgalactosaminyltransferase" );

is( $obj->definition( "Catalysis of ..." ), "Catalysis of ..." );
is( $obj->definition(), "Catalysis of ..." );

is( $obj->version( "666" ), "666" );
is( $obj->version(), "666" );

ok( $obj->ontology( "category 1 name" ) );
is( $obj->ontology()->name(), "category 1 name" );

my $ont = Bio::Ontology::Ontology->new();
ok( $ont->name( "category 2 name" ) );

ok( $obj->ontology( $ont ) );
is( $obj->ontology()->name(), "category 2 name" );

is( $obj->is_obsolete( 1 ), 1 );
is( $obj->is_obsolete(), 1 );

is( $obj->comment( "Consider the term ..." ), "Consider the term ..." );
is( $obj->comment(), "Consider the term ..." );

is( $obj->get_synonyms(), 0 );

$obj->add_synonym( ( "AA", "AB" ) );
my @al1 = $obj->get_synonyms();
is( scalar(@al1), 2 );
is( $al1[ 0 ], "AA" );
is( $al1[ 1 ], "AB" );

my @al2 = $obj->remove_synonyms();
is( $al2[ 0 ], "AA" );
is( $al2[ 1 ], "AB" );

is( $obj->get_synonyms(), 0 );
is( $obj->remove_synonyms(), 0 );

$obj->add_synonyms( ( "AA", "AB" ) );

is( $obj->identifier(undef), undef );
is( $obj->name(undef), undef );
is( $obj->definition(undef), undef );
is( $obj->is_obsolete(0), 0 );
is( $obj->comment(undef), undef );


$obj = Bio::Ontology::Term->new( 
    -identifier  => "0016847",
    -name        => "1-aminocyclopropane-1-carboxylate synthase",
    -definition  => "Catalysis of ...",
    -is_obsolete => 0,
    -version     => "6.6.6",
    -ontology    => "cat",
    -comment     => "X",
    -dbxrefs    => [
        Bio::Annotation::DBLink->new(-database => 'db1'),
        Bio::Annotation::DBLink->new(-database => 'db2')
    ],
    -references => []
);  

is( $obj->identifier(), "0016847" );
is( $obj->name(), "1-aminocyclopropane-1-carboxylate synthase" );
is( $obj->definition(), "Catalysis of ..." );
is( $obj->is_obsolete(), 0);
is( $obj->comment(), "X" );
is( $obj->version(), "6.6.6" );
is( $obj->ontology()->name(), "cat" );
is( scalar($obj->get_dbxrefs), 2);
is( scalar($obj->get_references), 0);

# test object factory for terms
my $fact = Bio::Ontology::TermFactory->new();
$obj = $fact->create_object(-name => "some ontology term");
isa_ok $obj, "Bio::Ontology::TermI";
is ($obj->name, "some ontology term");

$fact->type("Bio::Ontology::GOterm");
$obj = $fact->create_object(-name => "some ontology term",
			    -identifier => "GO:987654");
isa_ok $obj, "Bio::Ontology::TermI";
isa_ok $obj, "Bio::Ontology::GOterm";
is ($obj->name, "some ontology term");
is ($obj->identifier, "GO:987654");

$fact->type("Bio::Annotation::OntologyTerm");
$obj = $fact->create_object(-name => "some ontology term",
			    -identifier => "GO:987654",
			    -ontology => "nonsense");
isa_ok $obj, "Bio::Ontology::TermI";
isa_ok $obj, "Bio::AnnotationI";
is ($obj->name, "some ontology term");
is ($obj->identifier, "GO:987654");
is ($obj->tagname, "nonsense");
