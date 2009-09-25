# -*-Perl-*- Test Harness script for Bioperl
# $Id: Sim4.t 11525 2007-06-27 10:16:38Z sendu $

use strict;
BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 27);
	
	use_ok('Bio::Tools::Sim4::Results');
}

my $sim4 = Bio::Tools::Sim4::Results->new(-file=> test_input_file('sim4.rev'), -estisfirst=>0);
ok ( $sim4, 'new Sim4 results instance') ;


my $exonset = $sim4->next_exonset;
my @exons = $exonset->sub_SeqFeature(); 

is (scalar(@exons), 10);

my $exon = 1;
is $exons[$exon]->est_hit()->seq_id(), 'HSHNCPA1';
like($exons[$exon]->seq_id(), qr/human/);
is $exons[$exon]->strand(), -1;
is $exons[$exon]->start(), 1048;
is $exons[$exon]->end(), 1117;
is $exons[$exon]->score, 93;
is $exons[$exon]->est_hit()->seqlength(), 1198;


$sim4 = Bio::Tools::Sim4::Results->new(-file=> test_input_file('sim4.for.for'), -estisfirst=>0);
ok ( $sim4, 'new Sim4 results instance') ;

$exonset = $sim4->next_exonset;
@exons = $exonset->sub_SeqFeature(); 

is (scalar(@exons), 4);

$exon = 1;
is $exons[$exon]->est_hit()->seq_id(), 'hs_est';
is $exons[$exon]->seq_id(), 'human';
is $exons[$exon]->strand(), 1;
is $exons[$exon]->start(), 1377;
is $exons[$exon]->end(), 1500;
is $exons[$exon]->score, 99;
is $exons[$exon]->est_hit()->seqlength(), 479;

ok($sim4->next_exonset);
@exons = $exonset->sub_SeqFeature();

is $exons[$exon]->est_hit()->seq_id(), 'hs_est';
is $exons[$exon]->seq_id(), 'human';
is $exons[$exon]->strand(), 1;
is $exons[$exon]->est_hit->start(), 120;
is $exons[$exon]->est_hit->end(), 243;
is $exons[$exon]->score, 99;
is $exons[$exon]->est_hit()->seqlength(), 479;
