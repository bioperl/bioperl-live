# This is -*-Perl-*- code
# $Id $
# Bioperl Interface Test Harness for the Bio::SeqFeature::TypeI Interface
# The environment variable IMPLEMENTATION_NAME should be set to the interface
# implementation to be tested.

use strict;
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

## Test that we have a relrange object.
ok( $implementation->isa( 'Bio::RelRangeI' ) );

use Bio::PrimarySeq;
my $foo =
 Bio::PrimarySeq->new(
   '-id'     => 'foo',
   '-length' => 100
 );
ok( $foo );

use Bio::RelRange;
my $bar =
  Bio::RelRange->new(
    '-seq_id' => $foo,
    '-start'  => 6,
    '-end'    => 26,
    '-strand' => '+'
  );
ok( $bar );

my $baz =
  $implementation->new(
    '-seq_id' => $bar,
    '-start'  => 1,
    '-end'    => 3,
    '-strand' => '+'
  );
ok( $baz );

ok( $baz->low(),  1 );
ok( $baz->high(),    3 );
ok( $baz->strand(), 1 );
ok( $baz->abs_low(),  6 );
ok( $baz->abs_high(),    8 );
ok( $baz->abs_strand(), 1 );

my $sha =
  $implementation->new(
    '-seq_id' => $bar,
    '-start'  => 17,
    '-end'    => 21,
    '-strand' => -1
  );
ok( $sha );

ok( $sha->low( 'stranded' ),  1 );
ok( $sha->high( 'stranded' ),    5 );
ok( $sha->strand(), -1 );

## TODO: REMOVE
#warn "relative: $sha";
#warn "relative, plus: ".$sha->toRelRangeString( 0, 'plus' );
#warn "absolute: ".$sha->toRelRangeString( 1, 'stranded' );
#warn "absolute, plus: ".$sha->toRelRangeString( 1, 'plus' );

# Set the position policy for the tests
$sha->position_policy( 'plus' );
ok( $sha->position_policy(), 'plus' );

ok( $sha->abs_high(),   26 );
ok( $sha->abs_low(),    22 );
ok( $sha->abs_strand(), -1 );

# Now we should be able to get the stranded values again by overriding
# the object's position policy for a single call.
ok( $sha->abs_high( 'stranded' ),   79 );
ok( $sha->abs_low( 'stranded' ),    75 );

# Set the position policy for the tests
$sha->position_policy( 'stranded' );
ok( $sha->position_policy(), 'stranded' );

ok( $sha->abs_high(),   79 );
ok( $sha->abs_low(),    75 );

ok( $sha->abs_high( 'plus' ),   26 );
ok( $sha->abs_low( 'plus' ),    22 );


ok( $bar->abs2rel( 6 ), 1 );
ok( $bar->abs2rel( 26 ), 21 );
ok( $bar->rel2abs( 6 ), 11 );
ok( $bar->rel2abs( 26 ), 31 );
ok( $baz->abs2rel( 6 ), 1 );
ok( $baz->abs2rel( 26 ), 21 );
ok( $baz->rel2abs( 6 ), 11 );
ok( $baz->rel2abs( 26 ), 31 );

# Set the position policy for the tests
$sha->position_policy( 'plus' );
ok( $sha->position_policy(), 'plus' );
ok( $sha->abs2rel( 6 ), -15 );
ok( $sha->abs2rel( 26 ), 5 );
ok( $sha->rel2abs( 6 ), 27 );
ok( $sha->rel2abs( 26 ), 47 );

# Set the position policy for the tests
$sha->position_policy( 'stranded' );
ok( $sha->position_policy(), 'stranded' );
ok( $sha->abs2rel( 6 ), -68 );
ok( $sha->abs2rel( 26 ), -48 );
ok( $sha->rel2abs( 6 ), 80 );
ok( $sha->rel2abs( 26 ), 100 );

## TODO: ERE I AM.  This is a work-in-progress..

# Start out by setting the seq_id to something good.
#$implementation->seq_id( 'foobar' );
#$implementation->start( 1 );
#$implementation->end( 100 );
#$implementation->strand( 1 );
#
#my $subrange = $implementation->new();
#$subrange->seq_id( $implementation );
#$subrange->start( 1 );
#$subrange->end( 10 );
#$subrange->strand( -1 );
#$subrange->orientation_policy( 'dependent' );
#$subrange->position_policy( 'plus' );
#
## Test that the plus position policy does the right thing for the - strand.
#ok( $subrange->abs_low(), 100 );
#ok( $subrange->abs_high(), 91 );
#
## Test that the plus position policy does the right thing for the + strand.
#$subrange->strand( 1 );
#ok( $subrange->abs_low(), 1 );
#ok( $subrange->abs_high(), 10 );
