# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 12);

    use_ok('Bio::SeqFeature::Computation');
}

my ($comp_obj1, $comp_obj2, @sft);

ok $comp_obj1 = Bio::SeqFeature::Computation->new(
    -start => 1,
    -end   => 10,
);

is $comp_obj1->computation_id(332), 332, 'computation id';

ok $comp_obj1->add_score_value('P', 33), 'score value';

ok $comp_obj2 = Bio::SeqFeature::Computation->new(
    -start => 2,
    -end   => 10,
);

ok $comp_obj1->add_SeqFeature($comp_obj2, 'exon');
ok @sft = $comp_obj1->get_all_SeqFeature_types();

is $sft[0], 'exon', 'sft[0] is exon';


ok $comp_obj1 = Bio::SeqFeature::Computation->new(
    -start => 10,
    -end => 100,
    -strand => -1,
    -primary => 'repeat',
    -program_name => 'GeneMark',
    -program_date => '12-5-2000',
    -program_version => 'x.y',
    -database_name => 'Arabidopsis',
    -database_date => '12-dec-2000',
    -computation_id => 2231,
    -score    => { no_score => 334 },
);

is $comp_obj1->computation_id, 2231, 'computation id';
ok $comp_obj1->add_score_value('P', 33);
is ( ($comp_obj1->each_score_value('no_score'))[0], '334', 'score value');
