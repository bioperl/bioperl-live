## Bioperl Test Harness Script for Modules
##


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
BEGIN { $| = 1; print "1..2\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::Root::RootI;

print "ok 1\n";
$loaded = 1;

package MyObject;
@ISA = qw (Bio::Root::RootI);

sub new {
    my ($class) = shift;
    my $self = {};

    bless $self, $class;
    return $self;
}

package main;

$obj = MyObject->new();


eval {
     $obj->throw('an exception');
};

if( $@ =~ /an exception/ ) {
     print "ok 2\n";
} else {
  print "not ok 2\n";
  print STDERR "Exception $@\n";
}

