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

  $NUMTESTS = 2;
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

## Test that we have a type object.
ok( $implementation->isa( 'Bio::SeqFeature::TypeI' ) );

