# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 26);
	
	use_ok('Bio::Tools::Spidey::Results');
}

ok my $spidey = Bio::Tools::Spidey::Results->new(-file=> test_input_file('spidey.noalignment'),
                                                 -verbose => -1);

$spidey->close();

my $exonset = $spidey->next_exonset;
ok(!defined($exonset));

$spidey = Bio::Tools::Spidey::Results->new(-file=> test_input_file('spidey.test1'));
$exonset = $spidey->next_exonset;
my @exons = $exonset->sub_SeqFeature(); 
is @exons, 6;

is $spidey->genomic_dna_length(), 145732769;
is $spidey->splicesites(), 4;
is $spidey->est_coverage(), 100;
is $spidey->overall_percentage_id(), 99.7;
is $spidey->missing_mrna_ends(), 'neither';

is $exonset->seq_id(), 'lcl|chr2';
is $exonset->start(), 36356457;
is $exonset->end(), 36375798;

my $exon = 0;
is $exons[$exon]->est_hit()->seq_id(), 'lcl|tmpseq_0';
is $exons[$exon]->start(), 36375691;
is $exons[$exon]->end(), 36375798;
is $exons[$exon]->strand(), -1;
is $exons[$exon]->est_hit()->start(), 1;
is $exons[$exon]->est_hit()->end(), 108;
is $exons[$exon]->donor(), 1;
is $exons[$exon]->acceptor(), 0;

$exon = 1;
is $exons[$exon]->start(), 36369345;
is $exons[$exon]->end(), 36369492;
is $exons[$exon]->est_hit()->start(), 109;
is $exons[$exon]->est_hit()->end(), 256;
is $exons[$exon]->donor(), 0;
is $exons[$exon]->acceptor(), 1;

$spidey->close();
