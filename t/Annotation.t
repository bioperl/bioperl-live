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
    plan tests => 26;
}

use Bio::Annotation::Collection;
use Bio::Annotation::DBLink;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::OntologyTerm;

ok(1);

my $link1 = new Bio::Annotation::DBLink(-database => 'TSC',
					-primary_id => 'TSC0000030'
					);
ok $link1->isa('Bio::AnnotationI');
ok $link1->database(), 'TSC';
ok $link1->primary_id(), 'TSC0000030';

my $ac = Bio::Annotation::Collection->new();
ok $ac->isa('Bio::AnnotationCollectionI');

$ac->add_Annotation('dblink',$link1);
$ac->add_Annotation('dblink',
		    Bio::Annotation::DBLink->new(-database => 'TSC',
						 -primary_id => 'HUM_FABV'));

my $comment = Bio::Annotation::Comment->new( '-text' => 'sometext');
ok $comment->text, 'sometext';
$ac->add_Annotation('comment', $comment);


my $ref = Bio::Annotation::Reference->new( '-authors' => 'author line',
					   '-title'   => 'title line',
					   '-location'=> 'location line',
					   '-start'   => 12);
ok $ref->authors, 'author line';
ok $ref->title,  'title line';
ok $ref->location, 'location line';
ok $ref->start, 12;
ok $ref->database, 'MEDLINE';

$ac->add_Annotation('reference', $ref);


my $n = 0;
foreach my $link ( $ac->get_Annotations('dblink') ) {
    ok $link->database, 'TSC';
    ok $link->tagname(), 'dblink';
    $n++;
}
ok ($n, 2);

$n = 0;
my @keys = $ac->get_all_annotation_keys();
ok (scalar(@keys), 3);
foreach my $ann ( $ac->get_Annotations() ) {
    shift(@keys) if ($n > 0) && ($ann->tagname ne $keys[0]);
    ok $ann->tagname(), $keys[0];
    $n++;
}
ok ($n, 4);

$ac->add_Annotation($link1);

$n = 0;
foreach my $link ( $ac->get_Annotations('dblink') ) {
    ok $link->tagname(), 'dblink';
    $n++;
}
ok ($n, 3);

