# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# $Id$

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

    plan tests => 2 }
use Bio::Root::RootI;

ok(1);
package MyObject;
use vars qw(@ISA);
@ISA = qw (Bio::Root::RootI);

package main;

my $obj = MyObject->new();

eval {
     $obj->throw('an exception');
};

ok ( $@ =~ /an exception/);
