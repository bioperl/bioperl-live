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
    plan tests => 60;
}

use Bio::Annotation::Collection;
use Bio::Annotation::DBLink;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::OntologyTerm;
use Bio::Annotation::StructuredValue;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Cluster::UniGene;

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

# annotation of structured simple values (like swissprot's GN line)
my $ann = Bio::Annotation::StructuredValue->new();
ok ($ann->isa("Bio::AnnotationI"));

$ann->add_value([-1], "val1");
ok ($ann->value(), "val1");
$ann->value("compat test");
ok ($ann->value(), "compat test");
$ann->add_value([-1], "val2");
ok ($ann->value(-joins => [" AND "]), "compat test AND val2");
$ann->add_value([0], "val1");
ok ($ann->value(-joins => [" AND "]), "val1 AND val2");
$ann->add_value([-1,-1], "val3", "val4");
$ann->add_value([-1,-1], "val5", "val6");
$ann->add_value([-1,-1], "val7");
ok ($ann->value(-joins => [" AND "]), "val1 AND val2 AND (val3 AND val4) AND (val5 AND val6) AND val7");
ok ($ann->value(-joins => [" AND ", " OR "]), "val1 AND val2 AND (val3 OR val4) AND (val5 OR val6) AND val7");

$n = 1;
foreach ($ann->get_all_values()) {
    ok ($_, "val".$n++);
}

# nested collections
my $nested_ac = Bio::Annotation::Collection->new();
$nested_ac->add_Annotation('nested', $ac);

ok (scalar($nested_ac->get_Annotations()), 1);
($ac) = $nested_ac->get_Annotations();
ok $ac->isa("Bio::AnnotationCollectionI");
ok (scalar($nested_ac->get_all_Annotations()), 5);
$nested_ac->add_Annotation('gene names', $ann);
ok (scalar($nested_ac->get_Annotations()), 2);
ok (scalar($nested_ac->get_all_Annotations()), 6);
ok (scalar($nested_ac->get_Annotations('dblink')), 0);
my @anns = $nested_ac->get_Annotations('gene names');
ok $anns[0]->isa("Bio::Annotation::StructuredValue");
@anns = map { $_->get_Annotations('dblink');
	  } $nested_ac->get_Annotations('nested');
ok (scalar(@anns), 3);
ok (scalar($nested_ac->flatten_Annotations()), 2);
ok (scalar($nested_ac->get_Annotations()), 6);
ok (scalar($nested_ac->get_all_Annotations()), 6);

# OntologyTerm annotation
my $termann = Bio::Annotation::OntologyTerm->new(-label => 'test case',
						 -identifier => 'Ann:00001',
						 -category => 'dumpster');
ok ($termann->term);
ok ($termann->term->name, 'test case');
ok ($termann->term->identifier, 'Ann:00001');
ok ($termann->tagname, 'dumpster');
ok ($termann->category->name, 'dumpster');
ok ($termann->as_text, "dumpster|test case|Ann:00001");

# AnnotatableI
my $seq = Bio::Seq->new();
ok ($seq->isa("Bio::AnnotatableI"));
my $fea = Bio::SeqFeature::Generic->new();
ok ($fea->isa("Bio::AnnotatableI"));
my $clu = Bio::Cluster::UniGene->new();
ok ($clu->isa("Bio::AnnotatableI"));
