#!/usr/bin/perl -w

#----------------------------------------------------------------
# PROGRAM : filehandle.pl
# PURPOSE : To demonstrate passing filehandles between objects.
# AUTHOR  : Steve Chervitz (sac@bioperl.org)
# CREATED : 21 Mar 1997
# REVISION: $Id$
#
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#
# MODIFIED: sac --- Tue Mar 31 19:54:15 1998
#----------------------------------------------------------------

use lib "/home/steve/perl/lib";
use Bio::Root::Object  ();    
use Foo          ();

## In this demo, we want all output for this script to go to one file.
## To do this, we create a new output filehandle in the $GLOBAL object.
## Then, we set the output filehandles for all other objects created in this
## script to the global filehandle. 
##
## NOTE that this feature is considered experimental.
##  I haven't much occaision to use it. My gut feeling is that
##  having this feature adds unnecessary complexity, since I typically
##  need all output for a script to go to just one place.

# Create the GLOBAL object.
my $GLOBAL = new Bio::Root::Object(-NAME=>'global');

my $GOUT   = $GLOBAL->set_display(-WHERE=>'./filehandle.out');

### Uncomment the following line to redirect output to STDOUT.
#$GOUT = $GLOBAL->set_display(WHERE=>'-');

print $GOUT "\nObject Driver2.";
print $GOUT "\n-----------------\n";

$| = 1;

########################################
# Main

file_test();
 
print $GOUT "\nDone.\n\n";

exit;

#########################
sub file_test {

    print $GOUT "\n-------------------------------------------------------------\n";
    print $GOUT "Redirecting STDOUT to: ${\$GLOBAL->fh('name')}\n";
    print $GOUT "-------------------------------------------------------------\n\n";

    ### Create a new Foo object. 
    my $fooObj = new Foo(-NAME   =>'foo1', 
			 -FOO    =>150, 
			 -COLOR  =>'red', 
			 -BAR    =>250, 
			 -FLAVOR =>'lemon-lime' );

    $fooObj->display(-WHERE=>$GLOBAL->fh());
    ## Alternatively:  $fooObj->display(WHERE=>$GOUT); # since we defined $GOUT above.
    
    $fooObj->err() and $fooObj->print_err(-BEEP=>1);
}

#########################

