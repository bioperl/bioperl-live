# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);

eval {require Test::More;};
if ($@) {
	use lib 't/lib';
}
use Test::More;

eval {require Error;};
if ($@) {
	use lib 't/lib';
}
use_ok("Error");

use lib './examples/root/lib';

BEGIN {
	$NUMTESTS = 8;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	
    plan tests => $NUMTESTS; 
}

use Bio::Root::Exception;
use TestObject;
use Error qw(:try);

$Error::Debug = 1; 

# Set up a tester object.
ok my $test = TestObject->new();

is $test->data('Eeny meeny miney moe.'), 'Eeny meeny miney moe.';

# This demonstrates what will happen if a method defined in an
# interface that is not implemented in the implementating object.

eval { 
    try {
		$test->foo();
    }
    catch Bio::Root::NotImplemented with {
		my $err = shift;
		is ref $err, 'Bio::Root::NotImplemented';
    };
};

# TestObject::bar() deliberately throws a Bio::TestException, 
# which is defined in TestObject.pm
try {
    $test->bar;
}
catch Bio::TestException with {
    my $err = shift;
    is ref $err, 'Bio::TestException';
};


# Use the non-object-oriented syntax to throw a generic Bio::Root::Exception.
try {
    throw Bio::Root::Exception( "A generic error", 42 );
}
catch Bio::Root::Exception with {
    my $err = shift;
    is ref $err, 'Bio::Root::Exception';
    is $err->value, 42;
};

# Try to call a subroutine that doesn't exist. But because it occurs
# within a try block, the Error module will create a Error::Simple to
# capture it. Handy eh?

try {
	$test->foobar();
}
otherwise {
	my $err = shift;
	is ref $err, 'Error::Simple';
}; 

