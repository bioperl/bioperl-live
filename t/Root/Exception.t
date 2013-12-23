# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    eval {require Error;};

    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 8,
               -requires_module => 'Error');

    use lib './examples/root/lib';
    use_ok('TestObject');
}

use Error qw(:try);
$Error::Debug = test_debug();

# Set up a tester object.
ok my $test = TestObject->new(-verbose => test_debug());

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
