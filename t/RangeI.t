# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

## We start with some black magic to print on failure.
use vars qw(@funcs);
use Test;
use strict;
BEGIN {
  @funcs = qw(start end length strand);
  plan tests => 8;
}

use Bio::RangeI;

my $i = 1;
my $func;
while ($func = shift @funcs) {
  $i++;
  if(exists $Bio::RangeI::{$func}) {
    ok(1);
    eval {
      $Bio::RangeI::{$func}->();
    };
    ok( $@ );
  } else {
    ok(0);
  }
}
