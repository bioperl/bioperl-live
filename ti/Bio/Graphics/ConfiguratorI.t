# This is -*-Perl-*- code
# $Id $
# Bioperl Interface Test Harness for the Bio::Graphics::ConfiguratorI
# Interface. The environment variable IMPLEMENTATION_NAME should be set to the
# interface implementation to be tested.

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

  $NUMTESTS = 24;
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

## Test that we have a configurator object.
ok( $implementation->isa( 'Bio::Graphics::ConfiguratorI' ) );

## Test that there's no sections other than the default section to begin with.
ok( scalar( $implementation->get_sections() ), 0 );

## Test that we can set and then get a value in a new section.
$implementation->set( 'EST', 'fgcolor', 'chartreuse' );
ok( $implementation->get( 'EST', 'fgcolor' ), 'chartreuse' );

## Test that there's now one section.
ok( scalar( $implementation->get_sections() ), 1 );

## Test that we can set and then get a value in the default section.
$implementation->set( 'fgcolor', 'blue' );
ok( $implementation->get( 'fgcolor' ), 'blue' );

## Test that we can access it via the section 'general'
ok( $implementation->get( 'general', 'fgcolor' ), 'blue' );

## Test that we can access it via the section 'default'
ok( $implementation->get( 'default', 'fgcolor' ), 'blue' );

## Test that there's still only one section (the default one doesn't count).
ok( scalar( $implementation->get_sections() ), 1 );

## Test that we can set and then get a value in a another new section.
$implementation->set( 'dna', 'bgcolor', 'black' );
ok( $implementation->get( 'dna', 'bgcolor' ), 'black' );

## Test that there's now two sections.
ok( scalar( $implementation->get_sections() ), 2 );
## Test that they are in the right order.
my @sections = $implementation->get_sections();
ok( $sections[ 0 ], 'EST' );
ok( $sections[ 1 ], 'dna' );

## Test that we can set and then get another tag in an existing section.
$implementation->set( 'EST', 'bgcolor', 'mauve' );
ok( $implementation->get( 'EST', 'bgcolor' ), 'mauve' );

## Test that there's now two tags in that section.
ok( scalar( $implementation->get_tags( 'EST' ) ), 2 );

## Test that we can change a value and the set method returns the former value.
ok( $implementation->set( 'EST', 'bgcolor', 'pink' ), 'mauve' );
ok( $implementation->get( 'EST', 'bgcolor' ), 'pink' );

## Test that there's still two tags in that section.
ok( scalar( $implementation->get_tags( 'EST' ) ), 2 );

## Test that we can unset a tag by using the value 'undef'
$implementation->set( 'EST', 'bgcolor', 'undef' );
ok( !defined( $implementation->get( 'EST', 'bgcolor' ) ) );

## Test that there's now one tag in that section.
ok( scalar( $implementation->get_tags( 'EST' ) ), 1 );

## Test that we can remove a section by removing all of its tags.
$implementation->set( 'EST', 'fgcolor', 'undef' );
@sections = $implementation->get_sections();
ok( scalar( @sections ), 1 );
ok( $sections[ 0 ], 'dna' );

## Test setting and then getting a code segment (just a string, really).
$implementation->set( 'EST', 'height', 'sub { return "Short"; }' );
ok( $implementation->get( 'EST', 'height' ), 'sub { return "Short"; }' );

## Test get_and_eval
ok( $implementation->get_and_eval( 'EST', 'height' )->(), 'Short' );

#
###
### exercise the getters after setting
###
#
## exercise get_sections
#my @sections = $conf->get_sections();
#ok( $sections[0], 'EST' );
#
## exercise get_tags
#my @tags	= $conf->get_tags();
#ok( $tags[0], 'fgcolor' );
#
#@tags	= $conf->get_tags('EST');
#@tags	= sort(@tags);
#$tags	= join(" ", @tags);
#ok( $tags eq 'fgcolor height' );
#
## exercise get
#my @value	= $conf->get('fgcolor');
#ok( $value[0], 'blue' );
#@value	= $conf->get('EST','fgcolor');
#ok( $value[0], 'chartreuse' );
#@value	= $conf->get('EST','height');
#ok( $value[0], 'sub { return \"Short\"; }');
#
## exercise get_and_eval
#@value	= $conf->get('EST','height');
#ok( $value[0], 'Short' );

