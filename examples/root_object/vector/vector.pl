#!/usr/bin/perl -w

#---------------------------------------------------------------------
# PROGRAM : vector.pl
# PURPOSE : A basic driver script for testing Bio::Root::Vector.pm.
# AUTHOR  : Steve A. Chervitz
# CREATED : 4 May 1997 (sac@genome.stanford.edu)
# REVISION: $Id$
#
#  A PersonManager object is created which contains a Person data member.
#  The PersonManager need not have direct knowledge about
#  how many Persons it manages. All it needs to do is keep track of the most
#  recently added Person. The Person object keeps track of the rest through
#  the methods it inherits from Bio::Root::Vector.pm.
#
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#
# WARNING:
#  Vector.pm is considered experimental and is subject to change.
#  One outstanding problem is memory management, which is not handled
#  ideally in the current implementation. Scripts that manipulate many
#  Vector.pm objects may have a significant memory leak.
#---------------------------------------------------------------------

use lib "/home/steve/perl/lib";
use Bio::Root::Global qw(:devel);
use PersonManager     ();  # must be in the same dir as this script.

print "\nVector Driver.";
print "\n-----------------\n";

$| = 1;

########################################
# Main

#debug();

my $manager = new PersonManager(-name=>'Vector Tester');

print "\n\n-----------------------------------------------------\n";
print "   CREATING 10 Persons.";
print "\n-----------------------------------------------------\n\n";

my @names = ('joe','pat','peg','lou','al','bob','carol','chris','nat','ed');
for($i=0;$i<10;$i++) {
    $manager->add_person(-name=>$names[$i]);
}

$manager->display;

print STDERR "\nHit <RETURN> to continue..."; <STDIN>;

print "\n\n-----------------------------------------------------\n";
print "   SORTING BY NAME.";
print "\n-----------------------------------------------------\n\n";

$manager->sort_data('name');    
$manager->display;

print STDERR "\nHit <RETURN> to continue..."; <STDIN>;

print "\n\n-----------------------------------------------------\n";
print "   DELETING 'bob'.";
print "\n-----------------------------------------------------\n\n";

#debug(2);
$manager->remove_person('bob');

$manager->display;

print STDERR "\nHit <RETURN> to continue..."; <STDIN>;

$manager->clear_err();
#debug(0);


print "\n\n-----------------------------------------------------\n";
print "   DELETING FIRST.";
print "\n-----------------------------------------------------\n\n";

$manager->remove_person('first');

$manager->display;

print STDERR "\nHit <RETURN> to continue..."; <STDIN>;

print "\n\n-----------------------------------------------------\n";
print "   DELETING LAST.";
print "\n-----------------------------------------------------\n\n";

$manager->remove_person('last');

$manager->display;

print "\n\n-----------------------------------------------------\n";
print "   INSERTING 'fred' AFTER 'joe'.";
print "\n-----------------------------------------------------\n\n";

$manager->insert_person(-NEW=>'fred', -AFTER=>'joe');

$manager->display;

print "\n\n-----------------------------------------------------\n";
print "   CLONING FIRST PERSON (don't tell the NIH!).";
print "\n-----------------------------------------------------\n\n";

my $person = $manager->get_person('first')->clone;

printf( "%12s: %-12s\n", "SIZE OF ${\$person->name}", $person->size);
$person->display;

print "\n\n-----------------------------------------------------\n";
print "   VECTOR AFTER CLONING.";
print "\n-----------------------------------------------------\n\n";
$manager->display;

print "\nDone.\n\n";

#debug(2);

exit;

#########################

