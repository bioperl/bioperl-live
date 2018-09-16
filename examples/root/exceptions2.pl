#!/usr/bin/env perl

# This shows how Error.pm-based objects can be thrown 
# by Bio::Root::Root::throw() when Error.pm is available.
# When Error.pm isn't available, Bio::Root::Root::throw() 
# works as usual.
#
# It also demonstrates what happens when you use an outer eval{}
# instead of a try{} to trap thrown Error.pm-based exceptions. 
# The behavior is the same as when Error.pm is not used.
# This is important for backward compatibility.
#
# Author: Steve Chervitz <sac@bioperl.org>
#

use strict;

use lib qw(lib/ ../../);

# Uncomment this line to force Bio::Root::Root::throw() to 
# not use Error.pm even if it's available.
# Some of the tests in this script will be skipped .
#BEGIN { $main::DONT_USE_ERROR = 1; }

use Bio::Root::Root;
#use Bio::Root::Exception;  # Not necessary since Bio::Root::Root uses it.
use Error qw(:try);

my $foo = Bio::Root::Root->new();

if (!$main::DONT_USE_ERROR) {
    try {
        # This is the new, fancier way to handle exceptions. 
        # You must have Error.pm to do this (tarball included in this dir).
        
        print "[1] Throwing Error within try block via call to Bio::Root::Root::throw()\n";
        $foo->throw( -class => 'Bio::Root::Exception',
                     -text  => "Oopsie!",
                     -value => "123" 
                   );
    }
    
    catch Bio::Root::Exception with {
        my $err = shift;
        print "[1] Caught Bio::Root::Exception:\n$err";

    }

    otherwise {
        my $err = shift;
        print "[1] Caught other Error: ", ref($err), "\n$err";
    };

    
    print "\n\n";
}

eval {

    # This example demonstrates the traditional method for throwing
    # an exception using Bio::Root::Root->throw('string').
    # Notice how an exception of type Bio::Root::Exception is created.

    print "[2] Calling Bio::Root::Root->throw('string') within an eval{}\n";
    $foo->throw("Error message string.");

};

if($@) {
    print "[2] Caught eval{}-based exception: ", ref($@), "\n$@";
}
else {
    print "[2] Nothing to catch.\n";
}



print "\n\n";

eval {

    # This example shows that calling Error::throw directly within
    # an eval{} doesn't lead to a true value in $@ if
    # the error lacks a value. 

    print "[3] Attempting to throw a valueless Error within an eval{} block\n    (this should fail to be caught by Error.pm v0.13 but is caught by v0.14 and greater).\n";

    if( $ENV{OSTYPE} =~ /cygwin/ ) {
        die "[3] This causes a segmentation fault with cygwin perl! Skipping.\n";
    }

    throw Error::Simple ("A simple error.");

};

if($@) {
    print "[3] Caught eval{}-based exception: ", ref($@), "\n$@\n";
}
else {
    print "[3] Nothing to catch.\n";
}


print "\n\n";

eval {

    # This example shows that calling Error::throw directly within
    # an eval{} *does* lead to a true value in $@ if the error 
    # contains a non-zero value. 

    print "[4] Attempting to throw a valued Error within an eval{} block.\n";

    throw Error::Simple ("A simple error.", 42);

};

if($@) {
    print "[4] Caught eval{}-based exception: ", ref($@), "\n$@\n";
}
else {
    print "[4] Nothing to catch.\n";
}

print "\n\n";

if (!$main::DONT_USE_ERROR) {
    eval {

        # This example shows what happens if we try to create a
        # Bio::Root::IOException (a subclass of Bio::Root::Exception)
        # with a zero value. Bio::Root::Exception::new() catches this
        # faux pas and substitutes a value that will register as true in if($@).

        print "[5] Attempting to throw a zero-valued Bio::Root::IOException\n    within an eval{} block.\n";

        throw Bio::Root::IOException ( -text =>"An error with zero value.",
                                   -value => 0);

    };

    if($@) {
        print "[5] Caught eval{}-based zero-valued exception: ", ref($@), "\n$@\n";
    }
    else {
        print "[5] Nothing to catch.\n";
    }
    print "\n\n";
}


eval {

    # If Error::throw is called *indirectly* within an eval{}
    # (i.e., by calling a method which then calls Error::throw),
    # $@ is defined and it consists of a reference to the Error.pm object.

    print "[6] Attempting to throw Error indirectly within an eval{} block \n    via Bio::Root::Root::throw()\n";

    $foo->throw( -class => 'Bio::Root::Exception',
                 -text  => "Oopsie!",
                 -value => "456"
                );

};

if($@) {
    print "[6] Caught eval{}-based exception: ", ref($@), "\n$@";
}
else {
    print "[6] Nothing to catch.\n";
}

print "Done.\n";


