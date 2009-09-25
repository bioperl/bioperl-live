# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use warnings;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 62,
               -requires_module => 'Graph::Directed');
	
    use_ok('Bio::Ontology::GOterm');
    use_ok('Bio::Ontology::Ontology');
	use_ok('Bio::Annotation::DBLink');
}

my $obj = Bio::Ontology::GOterm->new();

# note that all *dblink* methods are now deprecated...

isa_ok( $obj,"Bio::Ontology::GOterm");

$obj->init();

like( $obj->to_string(), qr'-- GO id:' );

is( $obj->GO_id( "GO:0003947" ), "GO:0003947" );
is( $obj->GO_id(), "GO:0003947" );


is( $obj->get_dbxrefs(), 0 );

my @dblinks = (Bio::Annotation::DBLink->new(-primary_id => 'dAA'),
				 Bio::Annotation::DBLink->new(-primary_id => 'dAB'));

$obj->add_dbxref( -dbxrefs => \@dblinks );
is( scalar($obj->get_dbxrefs()), 2 );
my @df1 = $obj->get_dbxrefs();
is( $df1[ 0 ]->display_text, "dAA" );
is( $df1[ 1 ]->display_text, "dAB" );
is( $obj->get_dbxrefs(), 2 );

my @df2 = $obj->remove_dbxrefs();
is( $df2[ 0 ]->display_text, "dAA" );
is( $df2[ 1 ]->display_text, "dAB" );

is( $obj->get_dbxrefs(), 0 );
is( $obj->remove_dbxrefs(), 0 );


is( $obj->get_secondary_GO_ids(), 0 );

$obj->add_secondary_GO_id( ( "GO:0000000", "GO:1234567" ) );
is( scalar($obj->get_secondary_GO_ids()), 2 );
my @si1 = $obj->get_secondary_GO_ids();
is( $si1[ 0 ], "GO:0000000" );
is( $si1[ 1 ], "GO:1234567" );
is( $obj->get_secondary_GO_ids(), 2 );

my @si2 = $obj->remove_secondary_GO_ids();
is( $si2[ 0 ], "GO:0000000" );
is( $si2[ 1 ], "GO:1234567" );

is( $obj->get_secondary_GO_ids(), 0 );
is( $obj->remove_secondary_GO_ids(), 0 );



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

$obj->add_synonym( ( "AA", "AB" ) );
$obj->add_dbxref( -dbxrefs => \@dblinks );
$obj->add_secondary_GO_id( ( "GO:1234567", "GO:1234567" ) );

$obj->init();
is( $obj->identifier(), undef ); # don't make up identifiers
is( $obj->name(), undef );
is( $obj->definition(), undef );
is( $obj->is_obsolete(), 0 );
is( $obj->comment(), undef );


$obj = Bio::Ontology::GOterm->new( -go_id       => "0016847",
                                   -name        => "1-aminocyclopropane-1-carboxylate synthase",
                                   -definition  => "Catalysis of ...",
                                   -is_obsolete => 0,
                                   -version     => "6.6.6",
                                   -ontology    => "cat",
                                   -comment     => "X" );  

is( $obj->identifier(), "GO:0016847" );
is( $obj->name(), "1-aminocyclopropane-1-carboxylate synthase" );
is( $obj->definition(), "Catalysis of ..." );
is( $obj->is_obsolete(), 0 );
is( $obj->comment(), "X" );
is( $obj->version(), "6.6.6" );
is( $obj->ontology()->name(), "cat" );
