#!/usr/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM : destroy.pl
# PURPOSE : A basic driver script for testing destruction of Bio::Root::Object.pm
#           references and parent-child relationships.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 3 Nov 1996
# REVISION: $Id$
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#---------------------------------------------------------------------------

use lib "/home/steve/perl/lib";
use Bio::Root::Object ();    
use Bio::Root::Global qw(:std);    
use Foo               ();
use Outer             ();

print "\nParent-child driver.";
print "\n---------------------\n";

select(STDOUT); $| = 1;
my @errs = ();
my (@objs, $foo, $bar);

########################################
# Main

&parent_child;

&create_obj(2000);
print "Verify memory usage and hit <RETURN> to continue.";<STDIN>;

&destroy;
print "Verify memory usage and hit <RETURN> to continue.";<STDIN>;

&create_obj(2000);
print "Verify memory usage and hit <RETURN> to continue.";<STDIN>;

&test;
print "Verify memory usage and hit <RETURN> to continue.";<STDIN>;

exit 0;


#-----------------
sub create_obj {
#-----------------
# Create Foo objects

    my $num = shift || 1;

    print "\n-------------------------------------------------------------\n";
    print "Create Objects\n";
    print "-------------------------------------------------------------\n\n";

    foreach(1..$num) {
       push @objs, new Foo(-NAME   =>"foo$_",
			   -FOO    =>150, 
			   -BAR    =>250, 
			   -COLOR  =>'crimson',
			   -FLAVOR =>'lemon-lime',
			   -VERBOSE=> 1,    
			   );
    }
    print "\nCreated $num Foo objects.\n";
}


#-----------------
sub parent_child {
#-----------------
# Test the _drop_child() method (considered a 'protected' method; end-users 
# don't need to call this method in the normal course of affairs).
# Also shows the use of the global debug() function;

    print "\n-------------------------------------------------------------\n";
    print "Parent-Child Tester\n";
    print "-------------------------------------------------------------\n\n";

#    debug(1);

    my $foo = new Foo(-NAME   =>"BigFoo",
		      -FOO    =>15000, 
		      -BAR    =>250, 
		      -COLOR  =>'crimson',
		      -FLAVOR =>'lemon-lime',
		      -VERBOSE=> 1);

    print "\nAttempting to get the Foo object to drop its Bar child...\n";
    $bar = $foo->bar;

    $bar->parent->_drop_child($bar);
    printf "%s is still alive and knows its parent (%s), but it is orphaned.\n", $bar->to_string, $bar->parent->to_string;
    debug(0);

    print "\nHere's proof: attempt to drop the Bar child again...\n";
    eval { $bar->parent->_drop_child($bar); };
    if($@) {print $@; }

    debug(1);

    printf "\nAdding %s back to %s as an array data member\n",  $bar->to_string, $bar->parent->to_string;
    $bar->parent->{'array_member'}->[0] = $bar;
    $bar->parent->_drop_child($bar);
    
    printf "\nAdding %s back to %s as a hash data member\n",  $bar->to_string, $bar->parent->to_string;
    $bar->parent->{'hash_member'}->{'bar'} = $bar;
    $bar->parent->_drop_child($bar);

    debug(0);
    
    print "\nTry to drop the Bar child again...\n";
    eval { $bar->parent->_drop_child($bar); };
    if($@) {print $@; }

    # Add Bar back to Foo for further processing.
    $foo->{'Bar'} = $bar;
}


#-----------------
sub destroy {
#-----------------
    print "\n-------------------------------------------------------------\n";
    print "Destroy\n";
    print "-------------------------------------------------------------\n\n";

#    debug(2);
	if( ref $bar) {
#	    print "\nDESTROYING BAR....\n";
	    $bar->destroy;
	    undef $bar;
	} 

    while($_ = shift @objs) {
	if (ref $_) {
#	    $_->display;  # uncomment this to check for creation/destruction of IOManagers
#	    printf "\nDESTROYING OBJ #%s:  $_....\n", scalar @objs;
	    $_->destroy if ref $_;
	    undef $_;
#	    <STDIN>;
	}
    }
}


#-----------------
sub test {
#-----------------
    print "\n-------------------------------------------------------------\n";
    print "Test\n";
    print "-------------------------------------------------------------\n\n";

    print "Is foo alive?\n";
    print "Foo: $foo\nBar: $bar\n";

    print "\nAssigning onto foo.\n";
    $foo = "123456";

    print "Foo: $foo\n";

}
