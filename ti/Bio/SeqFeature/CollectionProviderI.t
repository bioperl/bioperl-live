# This is -*-Perl-*- code
# $Id $
# Bioperl Interface Test Harness for the Bio::SeqFeature::CollectionProviderI
# Interface. The environment variable IMPLEMENTATION_NAME should be set to the
# interface implementation to be tested.

use strict;
use Bio::Range;
use Bio::SeqFeature::SimpleCollection;
use Bio::SeqFeature::SimpleType;
use vars qw( $NUMTESTS $IMPLEMENTATION_NAME );

BEGIN { 
  # to handle systems with no installed Test module
  # we include the t dir (where a copy of Test.pm is located)
  # as a fallback
  eval { require Test; };
  if( $@ ) {
    use lib 't';
  }
  use Test;

  $NUMTESTS = 53;
  plan tests => $NUMTESTS;
}

my $testnum;
my $verbose = $ENV{ 'BIOPERLDEBUG' } || 0;

$IMPLEMENTATION_NAME = $ENV{ 'IMPLEMENTATION_NAME' };

my $implementation;

## This is always Test #1.
if( $IMPLEMENTATION_NAME ) {
  eval "require $IMPLEMENTATION_NAME; \$implementation = $IMPLEMENTATION_NAME->new();";
  print STDERR $@ if $@;
  ok( $implementation && not $@ );
} else {
  ok( 0 );
}

## End of black magic.
##
## Insert additional test code below but remember to change the
## $NUMTESTS variable to reflect the total number of tests that will
## be run.  Note that the require test (above) must be included in $NUMTESTS.

## Test that we have a collection provider object.
ok( $implementation->isa( 'Bio::SeqFeature::CollectionProviderI' ) );

## Test that the collection provider starts out empty.
ok( $implementation->get_collection()->feature_count() == 0 );

my $feature_0 =
  new Bio::SeqFeature::Generic(
    -unique_id => 'test_feature_0',
    -start => 10,
    -end => 100,
    -strand => -1
  );

# Test that we can add one
my $collection = new Bio::SeqFeature::SimpleCollection( $feature_0 );
ok( $collection->feature_count() == 1 );
$implementation->insert_collection( $collection );
ok( $implementation->get_collection()->feature_count() == 1 );

# Test that if we try to add the same one, we get an exception.
undef $@;
eval {
  $implementation->insert_collection( $collection );
};
ok( $@ );
ok( $implementation->get_collection()->feature_count() == 1 );

# Test we can update that one, with no exception.
undef $@;
eval {
  $implementation->update_collection( $collection );
};
ok( !$@ );
ok( $implementation->get_collection()->feature_count() == 1 );

my $feature_1 =
  new Bio::SeqFeature::Generic(
    -unique_id => 'test_namespace:test_feature_1',
    -display_name => 'test_namespace:test_feature_1_name',
    -start => 20,
    -end => 90,
    -strand => 0
  );

my $feature_2 =
  new Bio::SeqFeature::Generic(
    -unique_id => 'test_feature_2',
    -display_name => 'test_feature_2_name',
    -start => 10,
    -end => 100,
    -strand => 1
  );

# Test that we can add two
$collection = new Bio::SeqFeature::SimpleCollection( $feature_1, $feature_2 );
$implementation->insert_collection( $collection );
ok( $implementation->get_collection()->feature_count() == 3 );

my $feature_3 =
  new Bio::SeqFeature::Generic(
    -display_name => 'test_feature_3_name',
    -start => 1,
    -end => 1000,
    -strand => 0
  );

# Test that if we try to add two, one of which is new, we get an exception.
$collection = new Bio::SeqFeature::SimpleCollection( $feature_1, $feature_3 );
undef $@;
eval {
  $implementation->insert_collection( $collection );
};
ok( $@ );
# The interface does not specify what the state of the collection
# provider should be now.  It may or may not contain $feature_3.
# Adjust so that the state is deterministic again.
if( $implementation->get_collection()->feature_count() == 3 ) {
  $collection = new Bio::SeqFeature::SimpleCollection( $feature_3 );
  $implementation->insert_collection( $collection );
}

## Test that the collection now has 4 features.
ok( $implementation->get_collection()->feature_count() == 4 );

## Test that the get_collection() method returns the correct features.
use vars qw( %feature_hash );
undef %feature_hash;
foreach my $feature ( $implementation->get_collection()->features() ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by a single unique_id (that does exist, as a unique_id)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection( -unique_id => 'test_feature_2' )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by a single unique_id (that does exist, as a unique_id with a namespace)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -unique_id => 'test_feature_1',
    -namespace => 'test_namespace'
                                 )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by a single unique_id (that does NOT exist)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection( -unique_id => 'test_feature_3' )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by a single unique_id (that does exist, as a name)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection( -unique_id => 'test_feature_3_name' )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by multiple unique_ids (that do exist, one with a namespace)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -unique_ids => [ 'test_feature_1', 'test_feature_2' ],
    -namespace => 'test_namespace'
                                 )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by multiple unique_ids (that do not exist)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -unique_ids => [ 'fake', 'test_feature_3' ]
                                 )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by multiple unique_ids (one of which exists)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -unique_ids => [ 'fake', 'test_feature_2' ]
                                 )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by a single name (that does exist, as a name)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection( -name => 'test_feature_3_name' )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by a single name (that does exist, as a name with a namespace)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -name => 'test_feature_1_name',
    -namespace => 'test_namespace'
                                 )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by a single name (that does NOT exist)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection( -name => 'test_feature_0_name' )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by a single name (that does exist, as a unique_id)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection( -name => 'test_feature_0' )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by multiple names (that do exist, one with a namespace)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -names => [ 'test_feature_1_name', 'test_feature_3_name' ],
    -namespace => 'test_namespace'
                                 )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by multiple names (that do not exist)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -names => [ 'fake', 'test_feature_0_name' ]
                                 )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by multiple names (one of which exists)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -names => [ 'fake', 'test_feature_3_name' ]
                                 )->features()
                    ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by a single attribute
undef %feature_hash;
$feature_0->add_tag_value( 'sweetness', 'quite sweet' );
$feature_1->add_tag_value( 'sweetness', 'rather sweet' );
$feature_2->add_tag_value( 'sweetness', 'wicked sweet' );
foreach my $feature (
  $implementation->get_collection(
    -attribute => { 'sweetness' => 'quite sweet' }
  )->features()
) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by multiple attributes
undef %feature_hash;
$feature_0->add_tag_value( 'sourness', 'soso sour' );
$feature_1->add_tag_value( 'sourness', 'soso sour' );
$feature_3->add_tag_value( 'sourness', "ain't sour at all" );
foreach my $feature (
  $implementation->get_collection(
    -attribute => { 'sourness' => 'soso sour',
                    'sweetness' => 'rather sweet' }
  )->features()
) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'overlaps' rangetype, 'ignore' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'overlaps',
    -strandmatch => 'ignore',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => 0 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'overlaps' rangetype, 'ignore' strandmatch, again.
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'overlaps',
    -strandmatch => 'ignore',
    -range => new Bio::Range( -start => 1, -end => 9, -strand => 0 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'overlaps' rangetype, 'strong' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'overlaps',
    -strandmatch => 'strong',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => 1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'overlaps' rangetype, 'strong' strandmatch, 0 strand
## Curiously, the requirement is that this combination never returns anything.
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'overlaps',
    -strandmatch => 'strong',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => 0 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'overlaps' rangetype, 'weak' strandmatch
## This should return all features with no strand and with the given strand
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'overlaps',
    -strandmatch => 'weak',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => -1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'overlaps' rangetype, 'weak' strandmatch, 0 strand
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'overlaps',
    -strandmatch => 'weak',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => 0 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'contains' rangetype, 'weak' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'contains',
    -strandmatch => 'ignore',
    -range => new Bio::Range( -start => 20, -end => 90, -strand => -1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'contains' rangetype, 'ignore' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'contains',
    -strandmatch => 'ignore',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => -1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'contains' rangetype, 'weak' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'contains',
    -strandmatch => 'weak',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => -1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'contains' rangetype, 'strong' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'contains',
    -strandmatch => 'strong',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => -1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'contained_in' rangetype, 'ignore' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'contained_in',
    -strandmatch => 'ignore',
    -range => new Bio::Range( -start => 20, -end => 90, -strand => 1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'contained_in' rangetype, 'ignore' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'contained_in',
    -strandmatch => 'ignore',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => 1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'contained_in' rangetype, 'weak' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'contained_in',
    -strandmatch => 'weak',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => 1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'contained_in' rangetype, 'strong' strandmatch
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -rangetype => 'contained_in',
    -strandmatch => 'strong',
    -range => new Bio::Range( -start => 10, -end => 100, -strand => 1 )
  )->features()
) {
  $feature_hash{ $feature }++;
}
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering on a single type.
undef %feature_hash;
my $exon_type =
  new Bio::SeqFeature::SimpleType(
    -name => 'exon',
    -unique_id => 'exon',
    -description => 'a type used by the test for CollectionI objects'
  );
$feature_0->type( $exon_type );
foreach my $feature (
  $implementation->get_collection(
    -type => $exon_type
  )->features()
) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering on a multiple types.
undef %feature_hash;
my $sine_type =
  new Bio::SeqFeature::SimpleType(
    -name => 'sine',
    -unique_id => 'sine',
    -description => 'another type used by the test for CollectionI objects'
  );
$feature_1->type( $sine_type );
foreach my $feature (
  $implementation->get_collection(
    -types => [ $exon_type, $sine_type ]
  )->features()
) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering on an ancestor type.
undef %feature_hash;
my $transposon_type =
  new Bio::SeqFeature::SimpleType(
    -name => 'transposon',
    -unique_id => 'transposon',
    -description => 'yet another type used by the test for CollectionI objects'
  );
$transposon_type->add_child( $exon_type );
foreach my $feature (
  $implementation->get_collection(
    -type => $transposon_type
  )->features()
) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering on a string (unique_id of a type)
undef %feature_hash;
foreach my $feature (
  $implementation->get_collection(
    -type => 'transposon'
  )->features()
) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

# Test that we can remove one
$collection = new Bio::SeqFeature::SimpleCollection( $feature_1 );
$implementation->remove_collection( $collection );
ok( $implementation->get_collection()->feature_count() == 3 );

# Test that if we try to remove the same one, we get an exception
undef $@;
eval {
  $implementation->remove_collection( $collection );
};
ok( $@ );
ok( $implementation->get_collection()->feature_count() == 3 );

# Test that we can remove two
$collection = new Bio::SeqFeature::SimpleCollection( $feature_2, $feature_3 );
$implementation->remove_collection( $collection );
ok( $implementation->get_collection()->feature_count() == 1 );

# Test that when we remove two, one of which is present, we get an exception
$collection = new Bio::SeqFeature::SimpleCollection( $feature_1, $feature_0 );
undef $@;
eval {
  $implementation->remove_collection( $collection );
};
ok( $@ );
# The interface does not specify what the state of the collection
# provider should be now.  It may or may not contain $feature_0.
# Adjust so that the state is deterministic again.
if( $implementation->get_collection()->feature_count() == 1 ) {
  $collection = new Bio::SeqFeature::SimpleCollection( $feature_0 );
  $implementation->remove_collection( $collection );
}

## Test that the collection now has no features.
ok( $implementation->get_collection()->feature_count() == 0 );
