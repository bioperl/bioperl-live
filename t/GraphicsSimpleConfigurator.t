# This is -*-Perl-*- code
## Bioperl Test Harness Script for the Bio::Graphics::SimpleConfigurator Module
##
# $Id $

use strict;
use vars qw( $NUMTESTS $INTERFACE_NAME $IMPLEMENTATION_NAME $test_name );

BEGIN {

  $INTERFACE_NAME = "Bio::Graphics::ConfiguratorI";
  $IMPLEMENTATION_NAME = "Bio::Graphics::SimpleConfigurator";

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

  $NUMTESTS = 1;
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
