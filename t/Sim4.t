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
    plan tests => 129;
}

use Bio::Tools::Sim4::Results;
use Bio::Root::IO;
use Bio::SearchIO;

ok(1);
my $sim4 = new Bio::Tools::Sim4::Results(-file=> Bio::Root::IO->catfile("t","data","sim4.rev"), -estisfirst=>0);
ok $sim4;


my $exonset = $sim4->next_exonset;
my @exons = $exonset->sub_SeqFeature(); 

ok @exons, 10;
my $exon = 1;
ok $exons[$exon]->est_hit()->seq_id(), 'HSHNCPA1';
ok($exons[$exon]->seq_id(), qr/human/);
ok $exons[$exon]->strand(), -1;
ok $exons[$exon]->start(), 1048;
ok $exons[$exon]->end(), 1117;
ok $exons[$exon]->score, 93;
ok $exons[$exon]->est_hit()->seqlength(), 1198;


$sim4 = new Bio::Tools::Sim4::Results(-file=> Bio::Root::IO->catfile("t","data","sim4.for.for"), -estisfirst=>0);
ok $sim4;

$exonset = $sim4->next_exonset;
@exons = $exonset->sub_SeqFeature(); 

ok @exons, 4;
$exon = 1;
ok $exons[$exon]->est_hit()->seq_id(), 'hs_est';
ok $exons[$exon]->seq_id(), 'human';
ok $exons[$exon]->strand(), 1;
ok $exons[$exon]->start(), 1377;
ok $exons[$exon]->end(), 1500;
ok $exons[$exon]->score, 99;
ok $exons[$exon]->est_hit()->seqlength(), 479;

ok($sim4->next_exonset);
@exons = $exonset->sub_SeqFeature();

ok $exons[$exon]->est_hit()->seq_id(), 'hs_est';
ok $exons[$exon]->seq_id(), 'human';
ok $exons[$exon]->strand(), 1;
ok $exons[$exon]->est_hit->start(), 120;
ok $exons[$exon]->est_hit->end(), 243;
ok $exons[$exon]->score, 99;
ok $exons[$exon]->est_hit()->seqlength(), 479;


# new SearchIO parser for Sim4

# parse align format 0
my $parser = new Bio::SearchIO(-format => 'sim4',
			       -file   => 
			       Bio::Root::IO->catfile(qw(t data crypto.sim4-0))
			       );
my $r = $parser->next_result;
ok($r->query_name, 'cn416');
ok($r->query_length, 630);

my $hit = $r->next_hit;
ok($hit->name, 'Contig147');
ok($hit->description, 'Contig147.fa');
ok($hit->length, 1086);

my $hsp = $hit->next_hsp;
ok($hsp->query->start, 36);
ok($hsp->query->end, 132);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 191);
ok($hsp->hit->end, 286);
ok($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 133);
ok($hsp->query->end, 191);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 343);
ok($hsp->hit->end, 401);
ok($hsp->hit->strand, 1);

# parse align format 3
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data crypto.sim4-3))
			    );
$r = $parser->next_result;
ok($r->query_name, 'cn416');
ok($r->query_length, 630);
$hit = $r->next_hit;
ok($hit->name, 'Contig147');
ok($hit->description, 'Contig147.fa');
ok($hit->length, 1086);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 36);
ok($hsp->query->end, 132);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 191);
ok($hsp->hit->end, 286);
ok($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 133);
ok($hsp->query->end, 191);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 343);
ok($hsp->hit->end, 401);
ok($hsp->hit->strand, 1);

# parse align format 4
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data crypto.sim4-4))
			    );
$r = $parser->next_result;
ok($r->query_name, 'cn416');
ok($r->query_length, 630);

$hit = $r->next_hit;
ok($hit->name, 'Contig147');
ok($hit->description, 'Contig147.fa');
ok($hit->length, 1086);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 36);
ok($hsp->query->end, 132);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 191);
ok($hsp->hit->end, 286);
ok($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 133);
ok($hsp->query->end, 191);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 343);
ok($hsp->hit->end, 401);
ok($hsp->hit->strand, 1);


# do the other sim4 files
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data sim4.rev))
			    );
$r = $parser->next_result;
ok($r->query_name, '/nfs/disk21/birney/prog/wise2/example/human.rev');
ok($r->query_length, 5368);
$hit = $r->next_hit;
ok($hit->name, 'HSHNCPA1');
ok($hit->description, 'temp.cdna');
ok($hit->length, 1198);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 486);
ok($hsp->query->end, 503);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 10);
ok($hsp->hit->end, 27);
ok($hsp->hit->strand, -1);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 1048);
ok($hsp->query->end, 1117);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 194);
ok($hsp->hit->end, 265);
ok($hsp->hit->strand, -1);

# do the other sim4 files fwd
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data sim4.for.for))
			    );
$r = $parser->next_result;
ok($r->query_name, 'human.genomic');
ok($r->query_length, 5368);
$hit = $r->next_hit;
ok($hit->name, 'hs_est');
ok($hit->description, 'est.for');
ok($hit->length, 479);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 695);
ok($hsp->query->end, 813);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 1);
ok($hsp->hit->end, 119);
ok($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 1377);
ok($hsp->query->end, 1500);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 120);
ok($hsp->hit->end, 243);
ok($hsp->hit->strand, 1);

# do the other sim4 files fwd rev
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data sim4.for.rev))
			    );
$r = $parser->next_result;
ok($r->query_name, 'human.genomic');
ok($r->query_length, 5368);
$hit = $r->next_hit;
ok($hit->name, 'REVCOMP');
ok($hit->description, 'hn_est.rev');
ok($hit->length, 479);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 695);
ok($hsp->query->end, 813);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 1);
ok($hsp->hit->end, 119);
ok($hsp->hit->strand, -1);

$hsp = $hit->next_hsp;
ok($hsp->query->start, 1377);
ok($hsp->query->end, 1500);
ok($hsp->query->strand, 1);
ok($hsp->hit->start, 120);
ok($hsp->hit->end, 243);
ok($hsp->hit->strand, -1);
