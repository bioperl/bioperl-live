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
    plan tests => 52;
}

use Bio::Ontology::Term;
use Bio::Ontology::TermFactory;
  
my $obj = Bio::Ontology::Term->new();

ok( $obj->isa( "Bio::Ontology::Term" ) );

$obj->init();



ok( $obj->identifier( "0003947" ), "0003947" );
ok( $obj->identifier(), "0003947" );

ok( $obj->name( "N-acetylgalactosaminyltransferase" ), "N-acetylgalactosaminyltransferase" );
ok( $obj->name(), "N-acetylgalactosaminyltransferase" );

ok( $obj->definition( "Catalysis of ..." ), "Catalysis of ..." );
ok( $obj->definition(), "Catalysis of ..." );

ok( $obj->version( "666" ), "666" );
ok( $obj->version(), "666" );

ok( $obj->category( "category 1 name" ) );
ok( $obj->category()->name(), "category 1 name" );

my $cat = Bio::Ontology::Term->new();
ok( $cat->name( "category 2 name" ) );

ok( $obj->category( $cat ) );
ok( $obj->category()->name(), "category 2 name" );

ok( $obj->is_obsolete( 1 ), 1 );
ok( $obj->is_obsolete(), 1 );

ok( $obj->comment( "Consider the term ..." ), "Consider the term ..." );
ok( $obj->comment(), "Consider the term ..." );

ok( $obj->each_synonym(), 0 );

ok( $obj->add_synonyms( ( "AA", "AB" ) ) );
ok( $obj->each_synonym(), 2 );
my @al1 = $obj->each_synonym();
ok( $al1[ 0 ], "AA" );
ok( $al1[ 1 ], "AB" );
ok( $obj->each_synonym(), 2 );

my @al2 = $obj->remove_synonyms();
ok( $al2[ 0 ], "AA" );
ok( $al2[ 1 ], "AB" );

ok( $obj->each_synonym(), 0 );
ok( $obj->remove_synonyms(), 0 );

ok( $obj->add_synonyms( ( "AA", "AB" ) ) );



$obj->init();
ok( $obj->identifier(), undef );
ok( $obj->name(), undef );
ok( $obj->definition(), undef );
ok( $obj->is_obsolete(), 0 );
ok( $obj->comment(), undef );


$obj = Bio::Ontology::Term->new( -identifier  => "0016847",
                                 -name        => "1-aminocyclopropane-1-carboxylate synthase",
                                 -definition  => "Catalysis of ...",
                                 -is_obsolete => 0,
                                 -version     => "6.6.6",
                                 -category    => "cat",
                                 -comment     => "X" );  

ok( $obj->identifier(), "0016847" );
ok( $obj->name(), "1-aminocyclopropane-1-carboxylate synthase" );
ok( $obj->definition(), "Catalysis of ..." );
ok( $obj->is_obsolete(), 0 );
ok( $obj->comment(), "X" );
ok( $obj->version(), "6.6.6" );
ok( $obj->category()->name(), "cat" );

# test object factory for terms
my $fact = Bio::Ontology::TermFactory->new();
$obj = $fact->create_object(-name => "some ontology term");
ok $obj->isa("Bio::Ontology::TermI");
ok ($obj->name, "some ontology term");

$fact->type("Bio::Ontology::GOterm");
$obj = $fact->create_object(-name => "some ontology term",
			    -identifier => "GO:987654");
ok $obj->isa("Bio::Ontology::TermI");
ok (ref($obj), "Bio::Ontology::GOterm");
ok ($obj->name, "some ontology term");
ok ($obj->identifier, "GO:987654");

$fact->type("Bio::Annotation::OntologyTerm");
$obj = $fact->create_object(-name => "some ontology term",
			    -identifier => "GO:987654",
			    -category => "nonsense");
ok $obj->isa("Bio::Ontology::TermI");
ok $obj->isa("Bio::AnnotationI");
ok ($obj->name, "some ontology term");
ok ($obj->identifier, "GO:987654");
ok ($obj->tagname, "nonsense");
