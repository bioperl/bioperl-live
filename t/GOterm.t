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
    plan tests => 66;
}

use Bio::Ontology::GOterm;
use Bio::Ontology::Ontology;
  
my $obj = Bio::Ontology::GOterm->new();

ok( $obj->isa( "Bio::Ontology::GOterm" ) );

$obj->init();

ok( $obj->to_string() );


ok( $obj->GO_id( "GO:0003947" ), "GO:0003947" );
ok( $obj->GO_id(), "GO:0003947" );


ok( $obj->get_dblinks(), 0 );

ok( $obj->add_dblink( ( "dAA", "dAB" ) ) );
ok( $obj->get_dblinks(), 2 );
my @df1 = $obj->get_dblinks();
ok( $df1[ 0 ], "dAA" );
ok( $df1[ 1 ], "dAB" );
ok( $obj->get_dblinks(), 2 );

my @df2 = $obj->remove_dblinks();
ok( $df2[ 0 ], "dAA" );
ok( $df2[ 1 ], "dAB" );

ok( $obj->get_dblinks(), 0 );
ok( $obj->remove_dblinks(), 0 );


ok( $obj->get_secondary_GO_ids(), 0 );

ok( $obj->add_secondary_GO_id( ( "GO:-------", "1234567" ) ) );
ok( $obj->get_secondary_GO_ids(), 2 );
my @si1 = $obj->get_secondary_GO_ids();
ok( $si1[ 0 ], "GO:-------" );
ok( $si1[ 1 ], "GO:1234567" );
ok( $obj->get_secondary_GO_ids(), 2 );

my @si2 = $obj->remove_secondary_GO_ids();
ok( $si2[ 0 ], "GO:-------" );
ok( $si2[ 1 ], "GO:1234567" );

ok( $obj->get_secondary_GO_ids(), 0 );
ok( $obj->remove_secondary_GO_ids(), 0 );



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

ok( $obj->add_synonym( ( "AA", "AB" ) ) );
ok( $obj->get_synonyms(), 2 );
my @al1 = $obj->get_synonyms();
ok( $al1[ 0 ], "AA" );
ok( $al1[ 1 ], "AB" );
ok( $obj->get_synonyms(), 2 );

my @al2 = $obj->remove_synonyms();
ok( $al2[ 0 ], "AA" );
ok( $al2[ 1 ], "AB" );

ok( $obj->get_synonyms(), 0 );
ok( $obj->remove_synonyms(), 0 );



ok( $obj->add_synonym( ( "AA", "AB" ) ) );
ok( $obj->add_dblink( ( "dAA", "dAB" ) ) );
ok( $obj->add_secondary_GO_id( ( "GO:1234567", "GO:1234567" ) ) );


$obj->init();
ok( $obj->identifier(), "GO:-------" );
ok( $obj->name(), undef );
ok( $obj->definition(), undef );
ok( $obj->is_obsolete(), 0 );
ok( $obj->comment(), undef );


$obj = Bio::Ontology::GOterm->new( -go_id       => "0016847",
                                   -name        => "1-aminocyclopropane-1-carboxylate synthase",
                                   -definition  => "Catalysis of ...",
                                   -is_obsolete => 0,
                                   -version     => "6.6.6",
                                   -ontology    => "cat",
                                   -comment     => "X" );  

ok( $obj->identifier(), "GO:0016847" );
ok( $obj->name(), "1-aminocyclopropane-1-carboxylate synthase" );
ok( $obj->definition(), "Catalysis of ..." );
ok( $obj->is_obsolete(), 0 );
ok( $obj->comment(), "X" );
ok( $obj->version(), "6.6.6" );
ok( $obj->ontology()->name(), "cat" );

