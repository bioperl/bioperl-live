# This is -*-Perl-*- code
## Bioperl Test Harness Script for the Bio::SeqFeature::SimpleType Module
##
# $Id $

use strict;
use vars qw( $NUMTESTS $INTERFACE_NAME $IMPLEMENTATION_NAME $test_name );

BEGIN {

  $INTERFACE_NAME = "Bio::SeqFeature::TypeI";
  $IMPLEMENTATION_NAME = "Bio::SeqFeature::SimpleType";

  use Test::Harness;
  $ENV{ 'IMPLEMENTATION_NAME' } = $IMPLEMENTATION_NAME;
  use lib '..';

  $test_name = $INTERFACE_NAME;
  $test_name =~ s|::|/|g;
  $test_name = "ti/${test_name}.t";

  # to handle systems with no installed Test module
  # we include the t dir (where a copy of Test.pm is located)
  # as a fallback
  eval { require Test; };
  if( $@ ) {
      use lib 't';
  }
  use Test;

  $NUMTESTS = 120;
  plan tests => $NUMTESTS;
}

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the $NUMTESTS variable in the BEGIN block to reflect the
## total number of tests that will be run. 

print STDERR "Running interface test \"$test_name\":\n";
eval( "runtests( \"$test_name\" )" );
print STDERR $@ if $@;
ok( !$@ );

use Bio::SeqFeature::SimpleType;

# Test no args constructor
my $simpletype =
  new Bio::SeqFeature::SimpleType();
ok( $simpletype );

# Test the many-args constructor with multiple synonyms
$simpletype =
  new Bio::SeqFeature::SimpleType(
    -name => 'test_type',
    -unique_id => 'id0',
    -version => '1.0',
    -category => 'fake',
    -comment => 'this is fake',
    -synonyms => [ '0', 'test_type' ],
    -description => 'a type used by the test for SimpleType.pm'
  );
ok( $simpletype );
ok( $simpletype->name() eq 'test_type' );
ok( $simpletype->description() eq 'a type used by the test for SimpleType.pm' );
ok( $simpletype->unique_id() eq 'id0' );
ok( $simpletype->version() eq '1.0' );
ok( $simpletype->category() eq 'fake' );
ok( $simpletype->comment() eq 'this is fake' );
my @synonyms = $simpletype->each_synonym();
ok( scalar( @synonyms ) == 2 );
ok( ( $synonyms[ 0 ] eq '0' ) &&
    ( $synonyms[ 1 ] eq 'test_type' ) );

# Test the constructor with a single synonym
my $dummy_simpletype =
  new Bio::SeqFeature::SimpleType(
    -name => 'dummy',
    -unique_id => 'id-1',
    -synonyms => 'a synonym'
  );
ok( $dummy_simpletype );
@synonyms = $dummy_simpletype->each_synonym();
ok( scalar( @synonyms ) == 1 );
ok( $synonyms[ 0 ] eq 'a synonym' );

# Test that the is_* methods correctly return false right now.
ok( !$simpletype->is_child( $dummy_simpletype ) );
ok( !$simpletype->is_parent( $dummy_simpletype ) );
ok( !$simpletype->is_descendent( $dummy_simpletype ) );
ok( !$simpletype->is_ancestor( $dummy_simpletype ) );

# Test that the parents, children, ancestors, & descendents methods
# correctly return undef now.
ok( !defined( $simpletype->children() ) );
ok( !defined( $simpletype->parents() ) );
ok( !defined( $simpletype->descendents() ) );
ok( !defined( $simpletype->ancestors() ) );

#
# Create a four nodes to be used in further testing
#
my $type1 = new Bio::SeqFeature::SimpleType(
    -name => 'test-a-1',
    -unique_id => 'id-type-1',
  );
my $type2 = new Bio::SeqFeature::SimpleType(
    -name => 'test-b-2',
    -unique_id => 'id-type-2',
  );
my $type3 = new Bio::SeqFeature::SimpleType(
    -name => 'test-c-3',
    -unique_id => 'id-type-3',
  );
my $type4 = new Bio::SeqFeature::SimpleType(
    -name => 'test-d-4',
    -unique_id => 'id-type-4',
  );
ok( $type1 && $type2 && $type3 && $type4 );

#
# Test that we can add a child to type1 and then
# retrieve it using the children() method.
#
$type1->add_child($type2);
my @children = $type1->children();
ok( scalar( @children ) == 1 );
ok( $children[ 0 ]->name eq 'test-b-2' );

#
# Add the remaining nodes to the graph
#
$type1->add_child($type3);
$type2->add_child($type4);
$type3->add_child($type4);
@children = $type1->children();
ok( scalar( @children ) == 2 );
my @names = ( $children[ 0 ]->name(),
                    $children[ 1 ]->name() );
ok ( join( "", sort( @names ) ) eq 'test-b-2test-c-3' );
@children = $type2->children();
ok( scalar( @children ) == 1 );
ok( $children[ 0 ]->name() eq 'test-d-4' );
@children = $type3->children();
ok( scalar( @children ) == 1 );
ok( $children[ 0 ]->name() eq 'test-d-4' );

#
# test the parents method
#
ok( !defined( $type1->parents() ) );
my @parents = $type2->parents();
ok( scalar( @parents ) == 1 );
ok( $parents[ 0 ]->name() eq 'test-a-1' );
my @parents = $type3->parents();
ok( scalar( @parents ) == 1 );
ok( $parents[ 0 ]->name() eq 'test-a-1' );
my @parents = $type4->parents();
ok( scalar( @parents ) == 2 );
@names = ( $parents[ 0 ]->name(),
           $parents[ 1 ]->name() );
ok ( join( "", sort( @names ) ) eq 'test-b-2test-c-3' );


#
# Check the descendents method
#
ok( !defined( $type4->descendents() ) );
my @descendents =  $type1->descendents();
ok( scalar( @descendents ) == 3 );
@names = ( $descendents[ 0 ]->name(),
           $descendents[ 1 ]->name(),
           $descendents[ 2 ]->name() );
ok( join( "", sort( @names ) ) eq 'test-b-2test-c-3test-d-4' );
@descendents = $type2->descendents();
ok( scalar( @descendents ) == 1 );
ok( $descendents[ 0 ]->name() eq 'test-d-4' );
@descendents = $type3->descendents();
ok( scalar( @descendents ) == 1 );
ok( $descendents[ 0 ]->name() eq 'test-d-4' );

#
# Check the ancestors method
#
ok( !defined( $type1->ancestors() ) );
my @ancestors =  $type4->ancestors();
ok( scalar( @ancestors ) == 3 );
@names = ( $ancestors[ 0 ]->name(),
           $ancestors[ 1 ]->name(),
           $ancestors[ 2 ]->name() );
ok( join( "", sort( @names ) ) eq 'test-a-1test-b-2test-c-3' );
@ancestors = $type2->ancestors();
ok( scalar( @ancestors ) == 1 );
ok( $ancestors[ 0 ]->name() eq 'test-a-1' );
@ancestors = $type3->ancestors();
ok( scalar( @ancestors ) == 1 );
ok( $ancestors[ 0 ]->name() eq 'test-a-1' );



# 
# Test is_parent
#
ok( !$type1->is_parent( $type1 ) );
ok( !$type1->is_parent( $type2 ) );
ok( !$type1->is_parent( $type3 ) );
ok( !$type1->is_parent( $type4 ) );

ok(  $type2->is_parent( $type1 ) );
ok( !$type2->is_parent( $type2 ) );
ok( !$type2->is_parent( $type3 ) );
ok( !$type2->is_parent( $type4 ) );

ok(  $type3->is_parent( $type1 ) );
ok( !$type3->is_parent( $type2 ) );
ok( !$type3->is_parent( $type3 ) );
ok( !$type3->is_parent( $type4 ) );

ok( !$type4->is_parent( $type1 ) );
ok(  $type4->is_parent( $type2 ) );
ok(  $type4->is_parent( $type3 ) );
ok( !$type4->is_parent( $type4 ) );

#
# Test is_child
#
ok( !$type1->is_child( $type1 ) );
ok(  $type1->is_child( $type2 ) );
ok(  $type1->is_child( $type3 ) );
ok( !$type1->is_child( $type4 ) );

ok( !$type2->is_child( $type1 ) );
ok( !$type2->is_child( $type2 ) );
ok( !$type2->is_child( $type3 ) );
ok(  $type2->is_child( $type4 ) );

ok( !$type3->is_child( $type1 ) );
ok( !$type3->is_child( $type2 ) );
ok( !$type3->is_child( $type3 ) );
ok(  $type3->is_child( $type4 ) );

ok( !$type4->is_child( $type1 ) );
ok( !$type4->is_child( $type2 ) );
ok( !$type4->is_child( $type3 ) );
ok( !$type4->is_child( $type4 ) );

#
# Test is_ancestor
#
ok( !$type1->is_ancestor( $type1 ) );
ok( !$type1->is_ancestor( $type2 ) );
ok( !$type1->is_ancestor( $type3 ) );
ok( !$type1->is_ancestor( $type4 ) );

ok(  $type2->is_ancestor( $type1 ) );
ok( !$type2->is_ancestor( $type2 ) );
ok( !$type2->is_ancestor( $type3 ) );
ok( !$type2->is_ancestor( $type4 ) );

ok(  $type3->is_ancestor( $type1 ) );
ok( !$type3->is_ancestor( $type2 ) );
ok( !$type3->is_ancestor( $type3 ) );
ok( !$type3->is_ancestor( $type4 ) );

ok(  $type4->is_ancestor( $type1 ) );
ok(  $type4->is_ancestor( $type2 ) );
ok(  $type4->is_ancestor( $type3 ) );
ok( !$type4->is_ancestor( $type4 ) );

#
# Test is_descendent
#
ok( !$type1->is_descendent( $type1 ) );
ok(  $type1->is_descendent( $type2 ) );
ok(  $type1->is_descendent( $type3 ) );
ok(  $type1->is_descendent( $type4 ) );

ok( !$type2->is_descendent( $type1 ) );
ok( !$type2->is_descendent( $type2 ) );
ok( !$type2->is_descendent( $type3 ) );
ok(  $type2->is_descendent( $type4 ) );

ok( !$type3->is_descendent( $type1 ) );
ok( !$type3->is_descendent( $type2 ) );
ok( !$type3->is_descendent( $type3 ) );
ok(  $type3->is_descendent( $type4 ) );

ok( !$type4->is_descendent( $type1 ) );
ok( !$type4->is_descendent( $type2 ) );
ok( !$type4->is_descendent( $type3 ) );
ok( !$type4->is_descendent( $type4 ) );

#
# Test using string for is_* matching
#
ok ( $type2->is_parent( 'id-type-1' ) );
ok ( $type2->is_child( 'id-type-4' ) );
ok ( !$type4->is_child( 'id-type-1' ) );
ok ( !$type2->is_child( 'not-defined' ) );


#

## TODO: ERE I AM
