# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($HAVEGRAPHDIRECTED $DEBUG $NUMTESTS);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;
	eval {require Graph::Directed;
			$HAVEGRAPHDIRECTED=1;
		};

	if ($@) {
		$HAVEGRAPHDIRECTED = 0;
		warn "Cannot run tests as Graph::Directed is not installed\n";
	}
	plan tests => ($NUMTESTS = 51);
}

END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Cannot complete Term tests',1);
	}
}

exit 0 unless $HAVEGRAPHDIRECTED;

require Bio::Ontology::Term;
use Bio::Ontology::TermFactory;
use Bio::Annotation::DBLink;
use Bio::Annotation::Reference;

my $obj = Bio::Ontology::Term->new();

ok( $obj->isa( "Bio::Ontology::TermI" ) );

ok( $obj->identifier( "0003947" ), "0003947" );
ok( $obj->identifier(), "0003947" );

ok( $obj->name( "N-acetylgalactosaminyltransferase" ), "N-acetylgalactosaminyltransferase" );
ok( $obj->name(), "N-acetylgalactosaminyltransferase" );

ok( $obj->definition( "Catalysis of ..." ), "Catalysis of ..." );
ok( $obj->definition(), "Catalysis of ..." );

ok( $obj->version( "666" ), "666" );
ok( $obj->version(), "666" );

ok( $obj->ontology( "category 1 name" ) );
ok( $obj->ontology()->name(), "category 1 name" );

my $ont = Bio::Ontology::Ontology->new();
ok( $ont->name( "category 2 name" ) );

ok( $obj->ontology( $ont ) );
ok( $obj->ontology()->name(), "category 2 name" );

ok( $obj->is_obsolete( 1 ), 1 );
ok( $obj->is_obsolete(), 1 );

ok( $obj->comment( "Consider the term ..." ), "Consider the term ..." );
ok( $obj->comment(), "Consider the term ..." );

ok( $obj->get_synonyms(), 0 );

$obj->add_synonym( ( "AA", "AB" ) );
my @al1 = $obj->get_synonyms();
ok( scalar(@al1), 2 );
ok( $al1[ 0 ], "AA" );
ok( $al1[ 1 ], "AB" );

my @al2 = $obj->remove_synonyms();
ok( $al2[ 0 ], "AA" );
ok( $al2[ 1 ], "AB" );

ok( $obj->get_synonyms(), 0 );
ok( $obj->remove_synonyms(), 0 );

$obj->add_synonyms( ( "AA", "AB" ) );

ok( $obj->identifier(undef), undef );
ok( $obj->name(undef), undef );
ok( $obj->definition(undef), undef );
ok( $obj->is_obsolete(0), 0 );
ok( $obj->comment(undef), undef );


$obj = Bio::Ontology::Term->new( 
    -identifier  => "0016847",
    -name        => "1-aminocyclopropane-1-carboxylate synthase",
    -definition  => "Catalysis of ...",
    -is_obsolete => 0,
    -version     => "6.6.6",
    -ontology    => "cat",
    -comment     => "X",
    -dblinks    => [
        Bio::Annotation::DBLink->new(-database => 'db1'),
        Bio::Annotation::DBLink->new(-database => 'db2')
    ],
    -references => []
);  

ok( $obj->identifier(), "0016847" );
ok( $obj->name(), "1-aminocyclopropane-1-carboxylate synthase" );
ok( $obj->definition(), "Catalysis of ..." );
ok( $obj->is_obsolete(), 0);
ok( $obj->comment(), "X" );
ok( $obj->version(), "6.6.6" );
ok( $obj->ontology()->name(), "cat" );
ok( scalar($obj->get_dblinks), 2);
ok( scalar($obj->get_references), 0);

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
			    -ontology => "nonsense");
ok $obj->isa("Bio::Ontology::TermI");
ok $obj->isa("Bio::AnnotationI");
ok ($obj->name, "some ontology term");
ok ($obj->identifier, "GO:987654");
ok ($obj->tagname, "nonsense");
