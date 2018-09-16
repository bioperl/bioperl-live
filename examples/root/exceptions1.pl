#!/usr/bin/env perl

# A simple tester script for demonstrating how to throw and catch
# Error.pm objects. It also shows how to define new types of
# Error.pm-based objects. 
#
# It relies on the tester modules TestObject.pm and TestInterface.pm
# which you should also look at.
#
# Note that Bio::Root::NotImplemented is a subclass of Error.pm 
# and is defined in Bio::Root::Exception.pm
#
# This code requires Graham Barr's Error.pm module available from CPAN.
#
# Author: Steve Chervitz <sac@bioperl.org>
#

use strict;
use lib qw(lib/ ../../);
use Error qw(:try);
use TestObject;
use Getopt::Long;

# Command-line options:
my $eg = 0;        # which example to run (a number 1-4)
my $help = 0;      # print usage info

# $Error::Debug is set to true by default in Bio::Root::Interface.
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

(!$eg || $help) and die "Usage: $0 -eg 1|2|3|4 [-nodebug] [-h]\nOptions:\n$options";

print $Error::Debug ? "Try a -nodebug option to supress stack trace." : "Verbose stacktrace off.";
print "\n\n";

# Set up a tester object.
my $test = TestObject->new();
$test->data('Eeny meeny miney moe.');

try {

    test_notimplemented( $test ) if $eg == 1;

    test_custom_error( $test ) if $eg == 2;

    test_simple_error() if $eg == 3;

    # This subroutine doesn't even exist. But because it occurs within a try block,
    # the Error module will create a Error::Simple to capture it. Handy eh?
    if(  $eg == 4 ) {
	print "Test #4: Calling an undefined subroutine.\n";
	test_foobar();
    }

    # We shouldn't see this stuff.
    print "----\n";
    print "----\n";
    print "Some other code within the try block after the last throw...\n";
    print "----\n";
    print "----\n";
}

# Multiple catch blocks to handle different types of errors:

catch Bio::Root::NotImplemented with {
    my $error = shift;
    print "\nCaught a Bio::Root::NotImplemented.\n",
      "  file  : ", $error->file, "\n",
      "  line  : ", $error->line, "\n",
      "  text  : ", $error->text, "\n",
      "  value : ", $error->value, "\n",
      "  object: ", ref($error->object), "\n";

    print "\nstacktrace:\n", $error->stacktrace, "\n";

    print "\nstringify:\n$error\n";
    # The above line is equivalent to this:
    #print "\nstringify:\n", $error->stringify, "\n";
}

catch Bio::TestException with {
    # Since we know what type of error we're getting,
    # we can extract more information about the offending object
    # which is retrievable from the error object.
    my $error = shift;
    print "\nCaught a Bio::TestException.\n",
      "  file  : ", $error->file, "\n",
      "  line  : ", $error->line, "\n",
      "  text  : ", $error->text, "\n",
      "  value : ", $error->value, "\n",
      "  object: ", ref($error->object), "\n",
      "  data  : ", $error->object->data, "\n";

    print "\nstacktrace:\n", $error->stacktrace, "\n";
    print "\nstringify:\n", $error->stringify, "\n";

}

otherwise {
    # This is a catch-all handler for any type of error not handled above.
    my $error = shift;
    print "\nCaught an other type of error: ", ref($error), "\n",
      "  file  : ", $error->file, "\n",
      "  line  : ", $error->line, "\n",
      "  text  : ", $error->text, "\n",
      "  value : ", $error->value, "\n",
      "  object: ", ref($error->object), "\n";

#    print "\nstack_trace_dump:\n", $error->stack_trace_dump(), "\n";

    print "\nstacktrace:\n", $error->stacktrace, "\n";

    print "\nstringify:\n$error\n";

};  # This semicolon is essential.

print "\nDone $0\n";

sub test_notimplemented {

    my $test = shift;
    # This demonstrates what will happen if a method defined in an interface 
    # that is not implemented in the implementating object.

    print "Test #1: Inducing a Bio::Root::NotImplemented exception from TestObject\n";

    $test->foo();
}


sub test_custom_error {

    my $test = shift;

    # TestObject::bar() deliberately throws a Bio::TestException, 
    # which is defined in TestObject.pm

    print "Test #2: Throwing a Bio::TestException exception from TestObject\n";

    $test->bar;

}


sub test_simple_error {

    # Error::Simple comes with Error.pm and can have only a string and a value.

    print "Test #3: Throwing a Error::Simple object\n";

    throw Error::Simple( "A simple error", 42 );
}


