#!/usr/bin/env perl

# This shows that Error objects can be subclassed into more specialized types.
# Bio::Root::FileOpenException is a subclass of  Bio::Root::IOException. 
#
# We can write a generic handler to trap any type of IOException
# or we could handle FileOpenExceptions explicitly.
#
# To demo, run this script without any arguments, then try it with an argument
# that doesn't correspond to any file on your system (e.g., foobar). 
# Then try running with a valid file name.
#
# This requires Graham Barr's Error.pm module available from CPAN.
#
# Author: Steve Chervitz <sac@bioperl.org>
#

use strict;
use lib qw(lib/ ../../);
use Error qw(:try);
use Bio::Root::Exception;

try {
   print "Starting try block.\n";
   my $file = shift @ARGV || throw Bio::Root::IOException(-text=>"No file supplied.");

   open ( IN, $file) || throw Bio::Root::FileOpenException(-text=>"Can't open file \"$file\"", -value=> $!); 

   print "Opened file $file\n";

}
catch Bio::Root::IOException with {
    # This handler deals with IOException or any of its subclasses.
    # We could also write a handler with a `catch Bio::Root::FileOpenException'.
    # Such a handler would appear before this one.
    my $e = shift;
    print "Caught IOException:\n\n$e";
}
finally {
     close IN;
};

print "\nDone.\n";
