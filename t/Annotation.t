## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..10\n"; }
END {print "not ok 1\n" unless $loaded;}

use Bio::Annotation;
use Bio::Annotation::DBLink;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 



$link1 = new Bio::Annotation::DBLink(-database => 'TSC',
				     -primary_id => 'TSC0000030'
				     );
print  "ok 2\n";

if( $link1->database() eq 'TSC') { 
    print "ok 3\n";
} else {
    print "not ok 3\n";
}

if( $link1->primary_id() eq 'TSC0000030') { 
    print "ok 4\n";
} else {
    print "not ok 4\n";
}

my $a = Bio::Annotation->new ( '-description'  => 'something');
print  "ok 5\n";


$a->add_DBLink($link1);
print  "ok 6\n";

foreach $link ( $a->each_DBLink ) {
    $link->primary_id;
    $link->database;
}
print  "ok 7\n";


if( $a->description ne 'something' ) {
    print "not ok 8\n";
} else {
    print "ok 8\n";
}


$comment = Bio::Annotation::Comment->new( '-text' => 'sometext');

if( $comment->text ne 'sometext' ) {
    print "not ok 9\n";
} else {
    print "ok 9\n";
}


$ref = Bio::Annotation::Reference->new( '-authors' => 'author line',
					'-title'   => 'title line',
					'-location' => 'location line');

if( $ref->authors ne 'author line' || 
    $ref->title   ne 'title line' ||
    $ref->location ne 'location line' ) {
    print "not ok 10\n";
} else {
    print "ok 10\n";
}

