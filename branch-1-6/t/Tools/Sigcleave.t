# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18);
	
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Tools::Sigcleave');
}

#load n-terminus of MGR5_HUMAN as test seq
my $protein = "MVLLLILSVLLLKEDVRGSAQSSERRVVAHMPGDIIIGALFSVHHQPTVDKVHERKCGAVREQYGI";

ok my $seq= Bio::PrimarySeq->new(-seq => $protein);

ok my $sig = Bio::Tools::Sigcleave->new();
ok $sig->seq($seq);
ok my $sout = $sig->seq;
is $sout->seq, $protein;
is $sig->threshold, 3.5;
is $sig->threshold(5), 5;
is $sig->matrix, 'eucaryotic';
is $sig->matrix('procaryotic'), 'procaryotic';
is $sig->matrix('eucaryotic'), 'eucaryotic';

like $sig->pretty_print, qr/Maximum score 7/;
ok my %results = $sig->signals;

is $results{9}, 5.2, "unable to get raw sigcleave results";


$sig = Bio::Tools::Sigcleave->new(-seq=>$protein,
				 -threshold=>5);
ok %results = $sig->signals;
is $results{9}, 5.2, "unable to get raw sigcleave results";
is $sig->result_count, 5;
