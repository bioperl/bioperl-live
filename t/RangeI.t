## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
use vars qw(@funcs);
BEGIN {
  @funcs = qw(start end length strand);
  $| = 1; print "1..", scalar(@funcs) + 1, "\n"; 
  use vars qw($loaded);
}
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::RangeI;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $i = 1;
my $func;
while ($func = shift @funcs) {
  $i++;
  if(exists $Bio::RangeI::{$func}) {
    print "ok $i\n";
    eval {
      $Bio::RangeI::{$func}->();
    };
    print "not ok $i\n" unless $@;
  } else {
    print "not ok $i\n";
  }
}
