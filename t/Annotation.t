# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($HAVEGRAPHDIRECTED $DEBUG $NUMTESTS);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    eval {require Graph::Directed; 
	  $HAVEGRAPHDIRECTED=1;
	  require Bio::Annotation::OntologyTerm; };
    if ($@) {
	$HAVEGRAPHDIRECTED = 0;
    }
    plan tests => ($NUMTESTS = 89);
}

use Bio::Annotation::Collection;
use Bio::Annotation::DBLink;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::Target;
use Bio::Annotation::AnnotationFactory;

ok(1);

use Bio::Annotation::StructuredValue;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Cluster::UniGene;



# simple value

ok my $simple = Bio::Annotation::SimpleValue->new(
						  -tagname => 'colour',
						  -value   => '1'
						 ), ;
ok $simple, 1;
ok $simple->value, 1;
ok $simple->tagname, 'colour';

ok $simple->value(0), 0;
ok $simple->value, 0;
ok $simple, 0;

# link

ok my $link1 = new Bio::Annotation::DBLink(-database => 'TSC',
					-primary_id => 'TSC0000030'
					);
ok $link1->isa('Bio::AnnotationI');
ok $link1->database(), 'TSC';
ok $link1->primary_id(), 'TSC0000030';
ok $link1->as_text, 'Direct database link to TSC0000030 in database TSC';
my $ac = Bio::Annotation::Collection->new();
ok $ac->isa('Bio::AnnotationCollectionI');

$ac->add_Annotation('dblink',$link1);
$ac->add_Annotation('dblink',
		    Bio::Annotation::DBLink->new(-database => 'TSC',
						 -primary_id => 'HUM_FABV'));

my $comment = Bio::Annotation::Comment->new( '-text' => 'sometext');
ok $comment->text, 'sometext';
ok $comment->as_text, 'Comment: sometext';
$ac->add_Annotation('comment', $comment);



ok my $target = new Bio::Annotation::Target(-target_id  => 'F321966.1',
                                         -start      => 1,
                                         -end        => 200,
                                         -strand     => 1,
					 );
ok $ac->add_Annotation('target', $target);


my $ref = Bio::Annotation::Reference->new( '-authors' => 'author line',
					   '-title'   => 'title line',
					   '-location'=> 'location line',
					   '-start'   => 12);
ok $ref->authors, 'author line';
ok $ref->title,  'title line';
ok $ref->location, 'location line';
ok $ref->start, 12;
ok $ref->database, 'MEDLINE';
ok $ref->as_text, 'Reference: title line';
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
ok (scalar(@keys), 4);
foreach my $ann ( $ac->get_Annotations() ) {
    shift(@keys) if ($n > 0) && ($ann->tagname ne $keys[0]);
    ok $ann->tagname(), $keys[0];
    $n++;
}
ok ($n, 5);

$ac->add_Annotation($link1);

$n = 0;
foreach my $link ( $ac->get_Annotations('dblink') ) {
    ok $link->tagname(), 'dblink';
    $n++;
}
ok ($n, 3);

# annotation of structured simple values (like swissprot''is GN line)
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
ok (scalar($nested_ac->get_all_Annotations()), 6);
$nested_ac->add_Annotation('gene names', $ann);
ok (scalar($nested_ac->get_Annotations()), 2);
ok (scalar($nested_ac->get_all_Annotations()), 7);
ok (scalar($nested_ac->get_Annotations('dblink')), 0);
my @anns = $nested_ac->get_Annotations('gene names');
ok $anns[0]->isa("Bio::Annotation::StructuredValue");
@anns = map { $_->get_Annotations('dblink');
	  } $nested_ac->get_Annotations('nested');
ok (scalar(@anns), 3);
ok (scalar($nested_ac->flatten_Annotations()), 2);
ok (scalar($nested_ac->get_Annotations()), 7);
ok (scalar($nested_ac->get_all_Annotations()), 7);

if( $HAVEGRAPHDIRECTED ) {
# OntologyTerm annotation
    my $termann = Bio::Annotation::OntologyTerm->new(-label => 'test case',
						     -identifier => 'Ann:00001',
						     -ontology => 'dumpster');
    ok ($termann->term);
    ok ($termann->term->name, 'test case');
    ok ($termann->term->identifier, 'Ann:00001');
    ok ($termann->tagname, 'dumpster');
    ok ($termann->ontology->name, 'dumpster');
    ok ($termann->as_text, "dumpster|test case|");
} else { 
    for (1..6 ) { 
	skip('Graph::Directed not installed cannot test Bio::Annotation::OntologyTerm module',1);
    }
}

# AnnotatableI
my $seq = Bio::Seq->new();
ok ($seq->isa("Bio::AnnotatableI"));
my $fea = Bio::SeqFeature::Generic->new();
ok ($fea->isa("Bio::AnnotatableI"));
my $clu = Bio::Cluster::UniGene->new();
ok ($clu->isa("Bio::AnnotatableI"));




# tests for Bio::Annotation::AnnotationFactory

ok my $factory = new Bio::Annotation::AnnotationFactory();


# defaults to SimpleValue
ok $ann = $factory->create_object(-value => 'peroxisome',
                                  -tagname => 'cellular component');
ok ref $ann, 'Bio::Annotation::SimpleValue';

ok  $factory->type('Bio::Annotation::OntologyTerm');

ok $ann = $factory->create_object(-name => 'peroxisome',
                                  -tagname => 'cellular component');
ok ref $ann, 'Bio::Annotation::OntologyTerm';

#ok $ann = $factory->create_object(-text => 'this is a comment');
#ok ref $ann, 'Bio::Annotation::Comment';


ok $factory->type('Bio::Annotation::Comment');
ok $ann = $factory->create_object(-text => 'this is a comment');
ok ref $ann, 'Bio::Annotation::Comment';


# factory guessing the type: Comment
$factory = new Bio::Annotation::AnnotationFactory();
ok $ann = $factory->create_object(-text => 'this is a comment');
ok ref $ann, 'Bio::Annotation::Comment';

# factory guessing the type: Target
$factory = new Bio::Annotation::AnnotationFactory();
ok $ann = $factory->create_object(-target_id => 'F1234', -start => 1, -end => 10);
ok ref $ann, 'Bio::Annotation::Target';

# factory guessing the type: OntologyTerm
$factory = new Bio::Annotation::AnnotationFactory();
ok $ann = $factory->create_object(-name => 'peroxisome',
                                  -tagname => 'cellular component');
ok ref $ann, 'Bio::Annotation::OntologyTerm';
