# -*-Perl-*- Test Harness script for Bioperl
# $Id$

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18);
	
    use_ok('Bio::Annotation::Relation');
}

use Bio::Seq;
use Bio::SeqFeature::Generic;

my $seq = Bio::Seq->new(-seq=>'ACTGTGGCGTCAACT',
                        -desc=>'Sample Bio::Seq object',
                        -alphabet => 'dna',
                        -is_circular => 1
                       );

my $f1 = Bio::SeqFeature::Generic->new(-start => 1,
                                       -end   => 5,
                                       -strand  => 1,
                                       -primary_tag => 'region',
                                       -display_name => 'foo');

my $f2 = Bio::SeqFeature::Generic->new(-start => 7,
                                       -end   => 12,
                                       -strand  => 1,
                                       -primary_tag => 'region',
                                       -display_name => 'bar');

for my $f ($f1, $f2) {
    $seq->add_SeqFeature($f)
}

ok(my $relation = Bio::Annotation::Relation->new(-from => $f1,
                                              -to   => $f2,
                                              -class => 'Bio::SeqFeatureI',
                                              -type => 'binding_site'));

$seq->annotation->add_Annotation('binding_site', $relation);

my ($rel) = $seq->annotation->get_Annotations('binding_site');

is($rel->as_text, 'foo <-> bar(binding_site)');
is($rel->relation_class, 'Bio::SeqFeatureI');
is($rel->is_directed, 0);
isa_ok($rel->from, 'Bio::SeqFeatureI');
isa_ok($rel->to, 'Bio::SeqFeatureI');
is($rel->type, 'binding_site');

ok($relation = Bio::Annotation::Relation->new(-to   => $f2,
                                              -class => 'Bio::SeqFeatureI',
                                              -type => 'binding_site'));

ok($f1->annotation->add_Annotation('binding_site', $relation));

($rel) = $f1->annotation->get_Annotations('binding_site');
is($rel->as_text, 'self <-> bar(binding_site)');
is($rel->relation_class, 'Bio::SeqFeatureI');
is($rel->is_directed, 0);
ok(!defined($rel->from));
isa_ok($rel->to, 'Bio::SeqFeatureI');
is($rel->type, 'binding_site');

# is_directed
$rel = Bio::Annotation::Relation->new(-from => $f1,
                                        -to   => $f2,
                                        -class => 'Bio::SeqFeatureI',
                                        -type => 'binding_site',
                                        -is_directed  => 1);

is($rel->as_text, 'foo -> bar(binding_site)');
is($rel->is_directed, 1);

exit;
