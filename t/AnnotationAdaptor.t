# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN { 
    eval { require Test::More; };
    if( $@ ) { 
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 22;
	use_ok('Bio::SeqFeature::Generic');
	use_ok('Bio::SeqFeature::AnnotationAdaptor');
	use_ok('Bio::Annotation::DBLink');
	use_ok('Bio::Annotation::Comment');
	use_ok('Bio::Annotation::SimpleValue');
}


my $feat = Bio::SeqFeature::Generic->new();
$feat->add_tag_value("tag1", "value of tag1");
$feat->add_tag_value("tag1", "another value of tag1");
$feat->add_tag_value("tag2", "some value for a tag");

my $link1 = new Bio::Annotation::DBLink(-database => 'TSC',
                                        -primary_id => 'TSC0000030',
                                        #cannot use "tag2", already claimed by Bio::Annotation::SimpleValue by above $feat->add_tag_value() call.
                                        -tagname => "tag3"
                                       );
$feat->annotation->add_Annotation($link1);

my $anncoll = Bio::SeqFeature::AnnotationAdaptor->new(-feature => $feat);

is($anncoll->get_num_of_annotations(), 4);
is(scalar($anncoll->get_all_annotation_keys()), 3);

my @anns = $anncoll->get_Annotations("tag1");
my @vals = $feat->each_tag_value("tag1");

is (scalar(@anns), scalar(@vals));
for(my $i = 0; $i < @anns; $i++) {
  is ($anns[$i]->value(), $vals[$i]);
}

@anns = $anncoll->get_Annotations("tag3");
my @fanns = $feat->annotation->get_Annotations("tag3");
@vals = $feat->each_tag_value("tag3");

is (scalar(@fanns), 1);
is (scalar(@anns), 1);
is (scalar(@vals), 1);
is ($anns[0]->primary_id(), $vals[0]);

is ($anns[0]->primary_id(), $fanns[0]->primary_id());

my $comment = Bio::Annotation::Comment->new( '-text' => 'sometext');
$anncoll->add_Annotation('comment', $comment);

@fanns = $feat->annotation->get_Annotations("comment");
is (scalar(@fanns), 1);
is ($fanns[0]->text(), "sometext");

my $tagval = Bio::Annotation::SimpleValue->new(-value => "boring value",
					       -tagname => "tag2");
$anncoll->add_Annotation($tagval);

@anns = $anncoll->get_Annotations("tag3");
@fanns = $feat->annotation->get_Annotations("tag3");
@vals = $feat->each_tag_value("tag3");

is (scalar(@fanns), 1);
is (scalar(@anns), 1);
is (scalar(@vals), 1);
is ($anns[0]->primary_id(), $vals[0]);

is ($anncoll->get_num_of_annotations(), 6);


