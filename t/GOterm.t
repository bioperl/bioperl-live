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
    plan tests => 57;
}

use Bio::Ontology::GOterm;
  
my $obj = Bio::Ontology::GOterm->new();

ok( $obj->isa( "Bio::Ontology::GOterm" ) );

$obj->init();

ok( $obj->to_string() );


ok( $obj->GO_id( "GO:0003947" ), "GO:0003947" );
ok( $obj->GO_id(), "GO:0003947" );


ok( $obj->each_definition_reference(), 0 );

ok( $obj->add_definition_references( ( "dAA", "dAB" ) ) );
ok( $obj->each_definition_reference(), 2 );
my @df1 = $obj->each_definition_reference();
ok( $df1[ 0 ], "dAA" );
ok( $df1[ 1 ], "dAB" );
ok( $obj->each_definition_reference(), 2 );

my @df2 = $obj->remove_definition_references();
ok( $df2[ 0 ], "dAA" );
ok( $df2[ 1 ], "dAB" );

ok( $obj->each_definition_reference(), 0 );
ok( $obj->remove_definition_references(), 0 );


ok( $obj->each_secondary_GO_id(), 0 );

ok( $obj->add_secondary_GO_ids( ( "GO:-------", "1234567" ) ) );
ok( $obj->each_secondary_GO_id(), 2 );
my @si1 = $obj->each_secondary_GO_id();
ok( $si1[ 0 ], "GO:-------" );
ok( $si1[ 1 ], "GO:1234567" );
ok( $obj->each_secondary_GO_id(), 2 );

my @si2 = $obj->remove_secondary_GO_ids();
ok( $si2[ 0 ], "GO:-------" );
ok( $si2[ 1 ], "GO:1234567" );

ok( $obj->each_secondary_GO_id(), 0 );
ok( $obj->remove_secondary_GO_ids(), 0 );



ok( $obj->ontology_id( "0003947" ), "0003947" );
ok( $obj->ontology_id(), "0003947" );

ok( $obj->name( "N-acetylgalactosaminyltransferase" ), "N-acetylgalactosaminyltransferase" );
ok( $obj->name(), "N-acetylgalactosaminyltransferase" );

ok( $obj->definition( "Catalysis of ..." ), "Catalysis of ..." );
ok( $obj->definition(), "Catalysis of ..." );

ok( $obj->is_obsolete( 1 ), 1 );
ok( $obj->is_obsolete(), 1 );

ok( $obj->comment( "Consider the term ..." ), "Consider the term ..." );
ok( $obj->comment(), "Consider the term ..." );

ok( $obj->each_alias(), 0 );

ok( $obj->add_aliases( ( "AA", "AB" ) ) );
ok( $obj->each_alias(), 2 );
my @al1 = $obj->each_alias();
ok( $al1[ 0 ], "AA" );
ok( $al1[ 1 ], "AB" );
ok( $obj->each_alias(), 2 );

my @al2 = $obj->remove_aliases();
ok( $al2[ 0 ], "AA" );
ok( $al2[ 1 ], "AB" );

ok( $obj->each_alias(), 0 );
ok( $obj->remove_aliases(), 0 );


ok( $obj->add_aliases( ( "AA", "AB" ) ) );
ok( $obj->add_definition_references( ( "dAA", "dAB" ) ) );
ok( $obj->add_secondary_GO_ids( ( "GO:1234567", "GO:1234567" ) ) );


$obj->init();
ok( $obj->ontology_id(), "GO:-------" );
ok( $obj->name(), "" );
ok( $obj->definition(), "" );
ok( $obj->is_obsolete(), 0 );
ok( $obj->comment(), "" );


$obj = Bio::Ontology::GOterm->new( -go_id       => "0016847",
                                   -name        => "1-aminocyclopropane-1-carboxylate synthase",
                                   -definition  => "Catalysis of ...",
                                   -is_obsolete => 0,
                                   -comment     => "" );  

ok( $obj->ontology_id(), "GO:0016847" );
ok( $obj->name(), "1-aminocyclopropane-1-carboxylate synthase" );
ok( $obj->definition(), "Catalysis of ..." );
ok( $obj->is_obsolete(), 0 );
ok( $obj->comment(), "" );


