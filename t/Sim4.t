# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

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
    plan tests => 10;
}

use Bio::Tools::Sim4::Results;
use Bio::Root::IO;

ok(1);
my $sim4 = new Bio::Tools::Sim4::Results(-file=> Bio::Root::IO->catfile("t","sim4.rev"), -estisfirst=>0);
ok $sim4;


my $exonset = $sim4->next_exonset;
my @exons = $exonset->sub_SeqFeature(); 

ok @exons, 10;
my $exon = 1;
ok $exons[$exon]->est_hit()->seqname(), 'HSHNCPA1';
ok $exons[$exon]->seqname() =~ /human/;
ok $exons[$exon]->strand(), -1;
ok $exons[$exon]->start(), 1048;
ok $exons[$exon]->end(), 1117;
ok $exons[$exon]->score, 93;
ok $exons[$exon]->est_hit()->seqlength(), 1198;
