# -*-Perl-*-
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
BEGIN { $| = 1; print "1..11\n"; }
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


sub test ($$;$) {
    my($num, $true,$msg) = @_;
    $msg = '' if !defined $msg;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$link1 = new Bio::Annotation::DBLink(-database => 'TSC',
				     -primary_id => 'TSC0000030'
				     );
test 2, defined $link1 && ref($link1) =~ /Bio::Annotation::DBLink/;

test 3, ( $link1->database() eq 'TSC');

test 4, ( $link1->primary_id() eq 'TSC0000030');

my $a = Bio::Annotation->new ( '-description'  => 'something');
test 5, defined $a && ref($a) =~ /Bio::Annotation/;


$a->add_DBLink($link1);
test 6, defined $a;

foreach $link ( $a->each_DBLink ) {
    $link->primary_id;
    $link->database;
}
test 7, "passed link";


test 8, ( $a->description eq 'something' );


$comment = Bio::Annotation::Comment->new( '-text' => 'sometext');

test 9, ( $comment->text eq 'sometext' );


$ref = Bio::Annotation::Reference->new( '-authors' => 'author line',
					'-title'   => 'title line',
					'-location'=> 'location line',
					'-start'   => 12);

test 10, ( $ref->authors eq 'author line' && 
	   $ref->title   eq 'title line' &&
	   $ref->location eq 'location line' &&
	   $ref->start == 12);

test 11, $ref->database eq 'MEDLINE';
