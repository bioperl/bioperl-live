#!/usr/bin/perl -w

#----------------------------------------------------------------
# PROGRAM  : io.pl
# PURPOSE  : To demonstrate passing filehandles between objects.
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED  : 21 Mar 1997
# REVISION : $Id$
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#----------------------------------------------------------------

use lib "/home/steve/perl/lib";
use Bio::Root::Object  ();    
use Foo                ();

## In this demo, we want all output for this script to go to one file.
## To do this, we create a new output filehandle in the $main object.
## Then, we set the output filehandles for all other objects created in this
## script to the main filehandle. 
##
## NOTE that this feature is considered experimental.
##  I haven't much occaision to use it. My gut feeling is that
##  this adds unnecessary complexity, since typically
##  all output for a script is sent to just one place.

$| = 1;

########################################
# Main

# Create the main object.
my $main = new Bio::Root::Object(-NAME=>'main');
my $GOUT = $main->set_display(-WHERE=>'./io.out');

print $GOUT "\nI/O demo for Bio::Root::Objects.";
print $GOUT "\n-----------------------------------\n";

### Uncomment the following line to redirect output to STDOUT.
#$GOUT = $main->set_display(WHERE=>'-');

file_test();
 
print $GOUT "\nDone.\n\n";

exit;

#########################
sub file_test {

    print $GOUT "\n-------------------------------------------------------------\n";
    print $GOUT "Redirecting STDOUT to: ${\$main->fh('name')}\n";
    print $GOUT "-------------------------------------------------------------\n\n";

    ### Create a new Foo object. 
    my $fooObj = new Foo(-NAME   =>'foo1', 
			 -FOO    =>150, 
			 -COLOR  =>'red', 
			 -BAR    =>250, 
			 -FLAVOR =>'lemon-lime' );

    $fooObj->display(-WHERE=>$main->fh());
    ## Alternatively:  $fooObj->display(WHERE=>$GOUT); # since we defined $GOUT above.
    
    $fooObj->err() and $fooObj->print_err(-BEEP=>1);
}

#########################

