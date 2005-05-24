use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located
    # as a fallback
    eval { require Test; };
    if( $@ ) {
		  	use lib 't';
    }
    use Test;
    plan tests => 28;
}

use Bio::Tools::Spidey::Results;
ok(1);

my $spidey = new Bio::Tools::Spidey::Results(-file=> Bio::Root::IO->catfile("t", "data",
"spidey.noalignment"));
ok $spidey;

$spidey->close();
ok(1);

my $exonset = $spidey->next_exonset;
ok(!defined($exonset));

$spidey = new Bio::Tools::Spidey::Results(-file=> Bio::Root::IO->catfile("t", "data",
"spidey.test1"));
$exonset = $spidey->next_exonset;
my @exons = $exonset->sub_SeqFeature(); 
ok @exons, 6;

ok $spidey->genomic_dna_length(), 145732769;
ok $spidey->splicesites(), 4;
ok $spidey->est_coverage(), 100;
ok $spidey->overall_percentage_id(), 99.7;
ok $spidey->missing_mrna_ends(), 'neither';

ok $exonset->seq_id(), 'lcl|chr2';
ok $exonset->start(), 36356457;
ok $exonset->end(), 36375798;

my $exon = 0;
ok $exons[$exon]->est_hit()->seq_id(), 'lcl|tmpseq_0';
ok $exons[$exon]->start(), 36375691;
ok $exons[$exon]->end(), 36375798;
ok $exons[$exon]->strand(), -1;
ok $exons[$exon]->est_hit()->start(), 1;
ok $exons[$exon]->est_hit()->end(), 108;
ok $exons[$exon]->donor(), 1;
ok $exons[$exon]->acceptor(), 0;

$exon = 1;
ok $exons[$exon]->start(), 36369345;
ok $exons[$exon]->end(), 36369492;
ok $exons[$exon]->est_hit()->start(), 109;
ok $exons[$exon]->est_hit()->end(), 256;
ok $exons[$exon]->donor(), 0;
ok $exons[$exon]->acceptor(), 1;

$spidey->close();
ok(1);
