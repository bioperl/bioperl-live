# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    plan tests => 16;
}

use Bio::Annotation;
use Bio::Annotation::DBLink;

ok(1);

my $link1 = new Bio::Annotation::DBLink(-database => 'TSC',
				     -primary_id => 'TSC0000030'
				     );
ok defined $link1;
ok $link1->isa('Bio::Annotation::DBLink');

ok $link1->database(), 'TSC';

ok $link1->primary_id(), 'TSC0000030';

my $a = Bio::Annotation->new ( '-description'  => 'something');
ok defined $a;
ok $a->isa('Bio::Annotation');


$a->add_DBLink($link1);
ok defined $a;

foreach my $link ( $a->each_DBLink ) {
    $link->primary_id;
    $link->database;
}
ok ($a->each_DBLink, 1);


ok $a->description, 'something';


my $comment = Bio::Annotation::Comment->new( '-text' => 'sometext');

ok $comment->text, 'sometext';

my $ref = Bio::Annotation::Reference->new( '-authors' => 'author line',
					   '-title'   => 'title line',
					   '-location'=> 'location line',
					   '-start'   => 12);

ok $ref->authors, 'author line';
ok $ref->title,  'title line';
ok $ref->location, 'location line';
ok $ref->start, 12;

ok $ref->database, 'MEDLINE';
