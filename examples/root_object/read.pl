#!/usr/bin/perl -w

#----------------------------------------------------------------
# PROGRAM  : read.pl
# PURPOSE  : To demonstrate the read() method of Bio::Root::Object.pm
#            Uses three tester files read.test1, read.test2, read.test3
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED  : 20 Jul 1998
# REVISION : $Id$
# USAGE    : read.pl < filename
#
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#
# TODO: Demonstrate supplying a function reference.
#----------------------------------------------------------------

use lib "/home/steve/perl/lib";
use Bio::Root::Object  ();    
use Bio::Root::Global qw(:std);
use FileHandle;

monitor(1);

$file = './read.test';
open( IN, $file) || die "Can't open $file: $!\n";
$obj = Bio::Root::Object->new(-name => 'file handle typeglob demo');

#debug(1);
$data = $obj->read(-handle => \*IN);
#debug(0);

printf "\nSTRING DATA FROM %s:\n%s\n\n", $obj->name, $data;

print "===========================================================\n";

$file2 = './read.test2';
$fh = new FileHandle $file2, 'r';  # 'r' is optional.
$obj2 = Bio::Root::Object->new(-name => 'FileHandle object demo');
#debug(1);
@data = $obj->read(-handle => $fh);
#debug(0);

printf "\nLIST DATA FROM %s:\n%s\n\n", $obj2->name, join( ", ", @data );

print "===========================================================\n";

$file3 = './read.test3';
$obj = Bio::Root::Object->new(-name => 'file name + rec separator demo');

#debug(1);
@data = $obj->read(-file => $file3,
		   -rec_sep => ">");

printf "\nLIST DATA FROM %s:\n%s\n\n", $obj->name, join("\nRECORD: ", @data);

print "===========================================================\n";

# To test this demo, supply input to this script from STDIN: 
# E.g., $object3.pl < file
$obj = Bio::Root::Object->new(-name => 'STDIN demo');

#debug(1);
print "\nReading from STDIN..(type something then hit <RETURN> and ^D)\n";
$data = $obj->read();

printf "\nSTRING DATA FROM %s:\n%s\n\n", $obj->name, $data;

