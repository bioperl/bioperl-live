#!/usr/bin/env perl

# This shows how the examples work when Error.pm isn't installed.
# It also shows how to supress using Error.pm if it is installed
# and you don't want to use it for some reason.
#
# Here we use the eval{} style exception handling that's currently
# in vogue trapping Bioperl exceptions.
#
# Author: Steve Chervitz <sac@bioperl.org>
#


# Setting this variable simulates not having Error.pm installed.
BEGIN { $DONT_USE_ERROR = 1; }

use strict;
use lib qw(lib/ ../../);
use TestObject;
use Getopt::Long;

# Command-line options:
my $eg = 0;        # which example to run (a number 1-4)
my $help = 0;      # print usage info
$Error::Debug = 1; # enables verbose stack trace 

GetOptions( "debug!" => \$Error::Debug,
	    "eg=s"   => \$eg,    
	    "h"      => \$help   
	  ); 

my $options = << "OPTS";
      -eg  1|2|3|4   Run a particular example
      -nodebug       Deactivate verbose stacktrace
      -h             Print this usage
OPTS

(!$eg || $help) and die "Usage: $0 -eg 1|2|3|4|5 [-nodebug] [-h]\nOptions:\n$options";

# Set up a tester object.
my $test = TestObject->new();
$test->data('Eeny meeny miney moe.');

eval {

    test_notimplemented( $test ) if $eg == 1;

    test_custom_error( $test ) if $eg == 2;

    test_simple_error() if $eg == 3;

    # This subroutine doesn't even exist. But because it occurs within a try block,
    # the Error module will create a Error::Simple to capture it. Handy eh?
    if(  $eg == 4 ) {
	print "Test #4: Calling an undefined subroutine.\n";
	test_foobar();
    }

    # Throwing an exception the traditional bioperl way.
    if(  $eg == 5 ) {
	print "Test #5: Creating a Bio::Root::Root object and calling throw('string').\n";
        my $obj = Bio::Root::Root->new();
        $obj->throw("Throwing string from Bio::Root::Root object.");
    }

    # We shouldn't see this stuff.
    print "----\n";
    print "----\n";
    print "Some other code within the try block after the last throw...\n";
    print "----\n";
    print "----\n";
};

if($@) {
    my $error = shift;
    print "\nAn exception occurred:\n$@\n";
}
else {
    print "\nNo exception occurred\n";
}

print "\nDone $0\n";

sub test_notimplemented {

    my $test = shift;
    # This demonstrates what will happen if a method defined in an interface 
    # that is not implemented in the implementation.

    print "Test #1: Inducing a Bio::Root::NotImplemented exception from TestObject\n";

    $test->foo();
}    


sub test_custom_error {

    my $test = shift;

    # TestObject::bar() deliberately throws a Bio::Root::TestError, 
    # which is defined in TestObject.pm

    print "Test #2: Throwing a Bio::TestException exception from TestObject\n";

    $test->bar;  

}


sub test_simple_error {

    # This example won't work without Error.pm installed.
    # It shows how setting $DONT_USE_ERROR = 1 
    # really does simulate the absence of Error.pm.
    # The exception should report something like:
    # "Can't locate object method "throw" via package "Error::Simple"

    # Error::Simple comes with Error.pm and can have only a string and a value.

    print "Test #3: Throwing a Error::Simple object\n\n";

    print "This should fail to find object method 'throw' via package 'Error::Simple'\n";
    print "because Error.pm is not available.\n\n";

    throw Error::Simple( "A simple error", 42 );
}


