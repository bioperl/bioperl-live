# This is -*-Perl-*- code
# $Id $
# Bioperl Interface Test Harness for the Bio::SeqFeature::CollectionI Interface
# The environment variable IMPLEMENTATION_NAME should be set to the interface
# implementation to be tested.

use strict;
use Bio::Range;
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

  $NUMTESTS = 44;
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
## be run.  Note that the use test (above) must be included in $NUMTESTS.

## Test that we have a collection object.
ok( $implementation->isa( 'Bio::SeqFeature::CollectionI' ) );

## Test that the collection starts out empty.
ok( scalar( $implementation->features() ) == 0 );

my $feature_0 =
  new Bio::SeqFeature::Generic(
    -unique_id => 'test_feature_0',
    -start => 10,
    -end => 100,
    -strand => -1
  );

# Test that we can add one, and that add_features returns the number added.
ok( $implementation->add_features( $feature_0 ) == 1 );

# Test that if we try to add the same one, add_features returns 0.
ok( $implementation->add_features( $feature_0 ) == 0 );

my $feature_1 =
  new Bio::SeqFeature::Generic(
    -unique_id => 'test_feature_1',
    -start => 20,
    -end => 90,
    -strand => 0
  );

my $feature_2 =
  new Bio::SeqFeature::Generic(
    -unique_id => 'test_feature_2',
    -start => 10,
    -end => 100,
    -strand => 1
  );

# Test that we can add two, and that add_features returns the number added.
ok( $implementation->add_features( $feature_1, $feature_2 ) == 2 );

my $feature_3 =
  new Bio::SeqFeature::Generic(
    -unique_id => 'test_feature_3',
    -start => 1,
    -end => 1000,
    -strand => 0
  );

# Test that when we add two, one of which is new, add_features returns 1.
ok( $implementation->add_features( $feature_2, $feature_3 ) == 1 );

## Test that the collection now has 4 features.
ok( scalar( $implementation->features() ) == 4 );

## Test that the collection now has 4 features, using feature_count().
ok( scalar( $implementation->feature_count() ) == 4 );

## Test that the features() method returns the correct list.
use vars qw( %feature_hash );
undef %feature_hash;
foreach my $feature ( $implementation->features() ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test that the get_SeqFeatures() method returns the correct list.
undef %feature_hash;
foreach my $feature ( $implementation->get_SeqFeatures() ) {
  $feature_hash{ $feature }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test returning an iterator from the features() method.
undef %feature_hash;
my $feature_iterator = $implementation->features( -iterator => 1 );
ok( ref( $feature_iterator ) &&
    $feature_iterator->isa( "Bio::SeqFeature::IteratorI" ) );
while( $feature_iterator->has_more_features() ) {
  $feature_hash{ $feature_iterator->next_feature() }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test returning an iterator using the get_feature_stream() method.
undef %feature_hash;
$feature_iterator = $implementation->get_feature_stream();
ok( ref( $feature_iterator ) &&
    $feature_iterator->isa( "Bio::SeqFeature::IteratorI" ) );
while( $feature_iterator->has_more_features() ) {
  $feature_hash{ $feature_iterator->next_feature() }++;
}
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test passing a callback to the features() method.
sub _callback {
  my $feature = shift;
  my $collection = shift;

  unless( defined $collection ) {
    print STDERR "WARNING: the \$collection value (2nd argument to the _callback function) is undef.";
  }
  $feature_hash{ $feature }++;

  return 1;
} # callback( $feature, $collection )
undef %feature_hash;
# test that the features method returns true when given the callback option
ok( $implementation->features( -callback => \&_callback ) );
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by a single attribute
undef %feature_hash;
$feature_0->add_tag_value( 'sweetness', 'quite sweet' );
$feature_1->add_tag_value( 'sweetness', 'rather sweet' );
$feature_2->add_tag_value( 'sweetness', 'bitter sweet' );
$feature_iterator =
  $implementation->features(
    -iterator => 1,
    -attribute => { 'sweetness' => 'quite sweet' }
  );
while( $feature_iterator->has_more_features() ) {
  $feature_hash{ $feature_iterator->next_feature() }++;
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
$feature_iterator =
  $implementation->features(
    -iterator => 1,
    -attribute => { 'sourness' => 'soso sour',
                    'sweetness' => 'rather sweet' }
  );
while( $feature_iterator->has_more_features() ) {
  $feature_hash{ $feature_iterator->next_feature() }++;
}
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'overlaps' rangetype, 'ignore' strandmatch
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'overlaps',
  -strandmatch => 'ignore',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => 0 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'overlaps' rangetype, 'ignore' strandmatch,
## again, this time using the features_in_range() method.
undef %feature_hash;
$implementation->features_in_range(
  new Bio::Range( -start => 1, -end => 9, -strand => 0 ),
  -callback => \&_callback,
  -rangetype => 'overlaps',
  -strandmatch => 'ignore'
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'overlaps' rangetype, 'strong' strandmatch
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'overlaps',
  -strandmatch => 'strong',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => 1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'overlaps' rangetype, 'strong' strandmatch, 0 strand
## Curiously, the requirement is that this combination never returns anything.
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'overlaps',
  -strandmatch => 'strong',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => 0 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'overlaps' rangetype, 'weak' strandmatch
## This should return all features with no strand and with the given strand
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'overlaps',
  -strandmatch => 'weak',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => -1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'overlaps' rangetype, 'weak' strandmatch, 0
## strand, with the overlapping_features() method.
undef %feature_hash;
$implementation->overlapping_features(
  -callback => \&_callback,
  -strandmatch => 'weak',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => 0 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'contains' rangetype, 'weak' strandmatch, using
## the contained_features(..) method.
undef %feature_hash;
$implementation->contained_features(
  -callback => \&_callback,
  -strandmatch => 'ignore',
  -range => new Bio::Range( -start => 20, -end => 90, -strand => -1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'contains' rangetype, 'ignore' strandmatch
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'contains',
  -strandmatch => 'ignore',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => -1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'contains' rangetype, 'weak' strandmatch
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'contains',
  -strandmatch => 'weak',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => -1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'contains' rangetype, 'strong' strandmatch
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'contains',
  -strandmatch => 'strong',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => -1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering by 'contained_in' rangetype, 'ignore' strandmatch
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'contained_in',
  -strandmatch => 'ignore',
  -range => new Bio::Range( -start => 20, -end => 90, -strand => 1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } != 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'contained_in' rangetype, 'ignore' strandmatch,
## using the contained_in() method.
undef %feature_hash;
$implementation->contained_in(
  -callback => \&_callback,
  -strandmatch => 'ignore',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => 1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'contained_in' rangetype, 'weak' strandmatch
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'contained_in',
  -strandmatch => 'weak',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => 1 )
);
# test that it behaved as expected.
ok( ( $feature_hash{ $feature_0 } == 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } != 0 ) &&
    ( $feature_hash{ $feature_3 } != 0 ) );

## Test filtering by 'contained_in' rangetype, 'strong' strandmatch
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -rangetype => 'contained_in',
  -strandmatch => 'strong',
  -range => new Bio::Range( -start => 10, -end => 100, -strand => 1 )
);
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
$implementation->features(
  -callback => \&_callback,
  -type => $exon_type
);
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
$implementation->features(
  -callback => \&_callback,
  -types => [ $exon_type, $sine_type ]
);
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
$implementation->features(
  -callback => \&_callback,
  -type => $transposon_type
);
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

## Test filtering on a string (unique_id of a type)
undef %feature_hash;
$implementation->features(
  -callback => \&_callback,
  -type => 'transposon'
);
ok( ( $feature_hash{ $feature_0 } != 0 ) &&
    ( $feature_hash{ $feature_1 } == 0 ) &&
    ( $feature_hash{ $feature_2 } == 0 ) &&
    ( $feature_hash{ $feature_3 } == 0 ) );

# Test that we can remove one, and that remove_features returns the
# number removed.
ok( $implementation->remove_features( $feature_1 ) == 1 );

# Test that if we try to remove the same one, remove_features returns 0.
ok( $implementation->remove_features( $feature_1 ) == 0 );

## Test that the collection now has 3 features.
ok( scalar( $implementation->features() ) == 3 );

## Test that the collection now has 3 features, using feature_count().
ok( scalar( $implementation->feature_count() ) == 3 );

# Test that we can remove two, and that remove_features returns the
# number removed.
ok( $implementation->remove_features( $feature_2, $feature_3 ) == 2 );

# Test that when we remove two, one of which is present,
# remove_features returns 1.
ok( $implementation->remove_features( $feature_0, $feature_3 ) == 1 );

## Test that the collection is now empty
ok( scalar( $implementation->features() ) == 0 );
