#!/usr/bin/perl

# This shows how Error.pm-based objects can be thrown 
# by Bio::Root::RootI::throw() when Error.pm is available.
# When Error.pm isn't available, Bio::Root::Root::throw() 
# works as usual.
#
# It also demonstrates what happens when you use an outer eval{}
# instead of a try{} to trap thrown Error.pm-based exceptions. 
# The behavior is the same as when Error.pm is not used.
# This is important for backward compatibility.
#
# Author: Steve Chervitz <steve_chervitz@affymetrix.com>
#
# $Id$

use strict;

use lib "../../";

# Uncomment this line to force Bio::Root::Root::throw() to 
# not use Error.pm even if it's available.
#BEGIN { $main::DONT_USE_ERROR = 1; }

use Bio::Root::Root;
use Bio::Root::Exception;
use Error qw(:try);

my $foo = Bio::Root::Root->new();

try {
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


print "\n";

eval {

    # This example demonstrates that Bio::Root::Root->throw() 
    # won't use Error.pm when called with a string.

    print "[2] Throwing error within an eval{} and passing a string to Bio::Root::Root::throw()\n";
    $foo->throw("Error message string.");

};

if($@) {
    print "[2] Caught eval{}-based exception: ", ref($@), "\n$@";
}
else {
    print "[2] Nothing to catch.\n";
}



print "\n";

eval {

    # This example shows that calling Error::throw directly within
    # an eval{} doesn't lead to anything getting added to $@,
    # so don't do this. Use die() or croak() instead.

    print "[3] Attempting to throw Error directly within an eval{} block\n";

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


print "\n";

eval {

    # If Error::throw is called *indirectly* within an eval{}
    # (i.e., by calling a method which then calls Error::throw),
    # $@ is defined and it consists of a reference to the Error.pm object.

    print "[4] Attempting to throw Error indirectly within an eval{} block \nvia Bio::Root::Root::throw()\n";

    $foo->throw( -class => 'Bio::Root::Exception',
                 -text  => "Oopsie!",
                 -value => "456"
		);

};

if($@) {
    print "[4] Caught eval{}-based exception: ", ref($@), "\n$@";
}
else {
    print "[4] Nothing to catch.\n";
}

print "Done.\n";


