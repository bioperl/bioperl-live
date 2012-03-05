# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 19);

    use_ok('Bio::SeqFeature::Generic');
    use_ok('Bio::SeqFeature::FeaturePair');
}


my ($feat, $feat2, $pair);

ok $pair = Bio::SeqFeature::FeaturePair->new();

ok $feat = Bio::SeqFeature::Generic->new(
    -start => 40,
    -end => 80,
    -strand => 1,
    -primary => 'exon',
    -source => 'internal',
    -display_name => 'my exon feature',
    -tag => {
        silly => 20,
        new => 1
    }
);

ok $feat2 = Bio::SeqFeature::Generic->new(
    -start => 400,
    -end => 440,
    -strand => 1,
    -primary => 'other',
    -source => 'program_a',
    -phase => 1,
    -tag => {
        silly => 20,
        new => 1
    }
);

ok $pair->feature1($feat);
ok $pair->feature2($feat2);

is $pair->feature1, $feat, 'feature1 of pair stored';
is $pair->feature2, $feat2, 'feature2 of pair stored';
is $pair->start, 40, 'feature start';
is $pair->end, 80, 'feature end';
is $pair->primary_tag, 'exon', 'primary tag';
is $pair->source_tag, 'internal', 'source tag';
is $pair->hstart, 400, 'hstart';
is $pair->hend, 440, 'hend';
is $pair->hprimary_tag, 'other', 'hprimary tag';
is $pair->hsource_tag, 'program_a', 'hsource tag';

ok $pair->invert;
is $pair->end, 440, 'inverted end';
