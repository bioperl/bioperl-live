# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 23);
	
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

my $link1 = Bio::Annotation::DBLink->new(-database => 'TSC',
                                        -primary_id => 'TSC0000030',
                                        -tagname => "tag2"
                                       );
$feat->annotation->add_Annotation($link1);

my $anncoll = Bio::SeqFeature::AnnotationAdaptor->new(-feature => $feat);

is($anncoll->get_num_of_annotations(), 4);
is(scalar($anncoll->get_all_annotation_keys()), 2);

my @anns = $anncoll->get_Annotations("tag1");
my @vals = $feat->get_tag_values("tag1");

is (scalar(@anns), scalar(@vals));
for(my $i = 0; $i < @anns; $i++) {
  is ($anns[$i]->value(), $vals[$i]);
}

@anns = $anncoll->get_Annotations("tag2");
my @fanns = $feat->annotation->get_Annotations("tag2");
@vals = $feat->get_tag_values("tag2");

is (scalar(@fanns), 1);
is (scalar(@anns), 2);
is (scalar(@vals), 1);
is ($anns[0]->value(), $vals[0]);

is ($anns[1]->primary_id(), $fanns[0]->primary_id());

my $comment = Bio::Annotation::Comment->new( '-text' => 'sometext');
$anncoll->add_Annotation('comment', $comment);

@fanns = $feat->annotation->get_Annotations("comment");
is (scalar(@fanns), 1);
is ($fanns[0]->text(), "sometext");

my $tagval = Bio::Annotation::SimpleValue->new(-value => "boring value",
					       -tagname => "tag2");
$anncoll->add_Annotation($tagval);

@anns = $anncoll->get_Annotations("tag2");
@fanns = $feat->annotation->get_Annotations("tag2");
@vals = $feat->get_tag_values("tag2");

is (scalar(@fanns), 1);
is (scalar(@anns), 3);
is (scalar(@vals), 2);
is ($anns[0]->value(), $vals[0]);
is ($anns[1]->value(), $vals[1]);
is ($anns[2]->primary_id(), $fanns[0]->primary_id());
