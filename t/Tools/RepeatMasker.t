# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 28);

    use_ok('Bio::Tools::RepeatMasker');
}

my $inputfilename = test_input_file('repeatmasker.fa.out');
my $parser = Bio::Tools::RepeatMasker->new(-file => $inputfilename);

{
    my $rpt = $parser->next_result;
    is ($rpt->feature1->seq_id, "contig11600");
    is ($rpt->feature1->start, 1337);
    is ($rpt->feature1->end, 1407);
    is ($rpt->feature1->strand, 1);
    is ($rpt->feature1->primary_tag, "Simple_repeat");
    is ($rpt->feature1->source_tag, "RepeatMasker");
    is (scalar $rpt->feature1->get_tag_values('Target'), 3);
    is ($rpt->feature2->seq_id, "(TTAGGG)n");
    is ($rpt->feature2->start, 2);
    is ($rpt->feature2->end, 76);
    is ($rpt->feature1->primary_tag, "Simple_repeat");
    is ($rpt->feature1->source_tag, "RepeatMasker");
    is (scalar $rpt->feature2->get_tag_values('Target'), 3);
}

$parser->next_result for 1,2,3;

{
    my $rpt = $parser->next_result;
    is ($rpt->feature1->seq_id, "SL2.30ch10");
    is ($rpt->feature1->start, 38849);
    is ($rpt->feature1->end, 38940);
    is ($rpt->feature1->strand, -1);
    is ($rpt->feature1->primary_tag, "LTR");
    is ($rpt->feature1->source_tag, "RepeatMasker");
    is_deeply ( [ $rpt->feature1->get_tag_values('Target') ], ['LTR_PGSC0003DMS000000301_448',10681,10775] );
    is ($rpt->feature2->seq_id, "LTR_PGSC0003DMS000000301_448");
    is ($rpt->feature2->start, 10681);
    is ($rpt->feature2->end, 10775);
    is ($rpt->feature2->strand, -1);
    is ($rpt->feature1->primary_tag, "LTR");
    is ($rpt->feature1->source_tag, "RepeatMasker");
    is_deeply ([$rpt->feature2->get_tag_values('Target')],['SL2.30ch10',38849,38940]);
}

