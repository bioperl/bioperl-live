#!/usr/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM : error.pl
# PURPOSE : A basic driver script for testing Object.pm exception handling
#           using Bio::Root::Object.pm and Bio::Root::Err.pm.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 28 Oct 1996
# REVISION: $Id$
#
# USAGE   : err.pl >& err.out  # Traps actual exceptions as thrown
#                              # in addition to reports created by script.
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#
# MODIFIED:
#  sac --- Wed Dec  2 02:58:54 1998
#     * Changed name to err.pl
#  sac --- Tue Mar 31 19:59:12 1998
#---------------------------------------------------------------------------

use lib "/home/steve/perl/lib";
use Bio::Root::Object ();    
use Bio::Root::Global qw(:std);    
use Foo               ();
use Outer             ();

print "\nException handling demo for Bio::Root::Object.pm";
print "\n---------------------------------------------------\n";

$| = 1;
my @errs = ();

########################################
# Main

#verbosity(1);   # Try this for extra detail in the warnings/exceptions.
#verbosity(-1);   # Try this for reduced detail in the warnings/exceptions.
#strictness(2);   # Strict = 2 converts nonfatal warnings into fatal exceptions.
#strictness(-2);   # Strict = -2 converts fatal exceptions into nonfatal warnings.

nonfatal_error();  
fatal_error();  
contained_error(); 
throw_error();  

print "\nDone.\n\n";

exit 0;


#########################
sub nonfatal_error {

    print "\n-------------------------------------------------------------\n";
    print "Nonfatal Error Demo (\"Warning\")\n";
    print "-------------------------------------------------------------\n\n";

    ## Error: forgot to specify Foo color. This type of error is considered recoverable,
    ##        since the Foo object will simply select a default color.
    ##  NOTE: Parameter tags can be upper or lower case.
    $fooObj = new Foo(-NAME   =>'foo1',
		      -FOO    =>150, 
		      -BAR    =>250, 
		      -FLAVOR =>'lemon-lime',
#		      -RECORD_ERR=> 1,  # Try this to store the warnings in the object,
                                        # accessible via $object->print_err() or $object->err().
			  );

    $fooObj->display(-HEADER=>1);
    
    print "FULL ERROR REPORT:\n";
    print "--------------------\n";
    $fooObj->err() 
	? $fooObj->print_err() 
	: print "\nfooObj foo okay.\n";
    
}


#########################
sub fatal_error {

    print "\n-------------------------------------------------------------\n";
    print "Fatal Error Demo\n";
    print "-------------------------------------------------------------\n\n";

    ## Error: forgot to specify foo data. This type of error is 
    ##        considered unrecoverable.

    eval { $fooObj = new Foo(-NAME   =>'foo1', 
			      -BAR    =>250, 
			      -FLAVOR =>'lemon-lime') };

    if($@) {
	print "\n\nOOPS! Something really bad happened!\nCan't use Foo object. Sorry.\n\n$@\n";

	# Reconstruct the Bio::Root::Err.pm object from the $@ string:
	my $err = new Bio::Root::Err($@);
	
	print "\n\nReconstructing the Err.pm object from the \$\@ string:\n";
	print "------------------------------------------------------------\n";
	printf "MESSAGE AND NOTE: \n%s\n\n", $err->string('msgnote');
	printf "CONTEXT: %s\n", $err->context();
	printf "STACK STRING: %s\n", $err->stack();

	return;
    }


}

#########################
sub contained_error {
    print "\n-------------------------------------------------------------\n";
    print "Contained Error Demo\n";
    print "-------------------------------------------------------------\n\n";
    
    print "Setting record_err to true\n\n";
    record_err(1); # Attaches Err objects to the objects that generate them.

    ## Error in a contained object: bad data for Bar object.
    my $outer = new Outer(-NAME   =>'outer1', 
			  -FOODAT =>[-NAME  =>'foo2', 
				     -COLOR =>'mauve', 
				     -FLAVOR =>'cherry', 
				     -FOO   =>150, 
				     -BAR   =>'23.1']
			  );
    
    $outer->display(-HEADER=>1);

    print "\n\n---------------------------------\n";
    print "Testing err_string:";
    print "\n---------------------------------\n";

    print $outer->err_string();
    record_err(0);

}


#########################
sub throw_error {

    print "\n-------------------------------------------------------------\n";
    print "Throw Error Demo\n";
    print "-------------------------------------------------------------\n\n";
    
    my $catcher = new Bio::Root::Object(-NAME=>'Catcher');

    eval { 
	my $fooObj = new Foo(-NAME   =>'foo1', 
			      -FOO    =>150, 
			      -COLOR  =>'green', 
			      -BAR    =>250, 
			      -FLAVOR =>'lemon-lime');

	$fooObj->change_color('brown'); 
    };
    ## Catch the exception with the catcher object and re-throw it.
    if($@) { $catcher->throw($@); }
}


