# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 130;
}

use_ok('Bio::Tools::Sim4::Results');
use_ok('Bio::Root::IO');
use_ok('Bio::SearchIO');

my $sim4 = new Bio::Tools::Sim4::Results(-file=> Bio::Root::IO->catfile("t","data","sim4.rev"), -estisfirst=>0);
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


$sim4 = new Bio::Tools::Sim4::Results(-file=> Bio::Root::IO->catfile("t","data","sim4.for.for"), -estisfirst=>0);
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


# new SearchIO parser for Sim4

# parse align format 0
my $parser = new Bio::SearchIO(-format => 'sim4',
			       -file   => 
			       Bio::Root::IO->catfile(qw(t data crypto.sim4-0))
			       );
my $r = $parser->next_result;
is ($r->query_name, 'cn416');
is ($r->query_length, 630);

my $hit = $r->next_hit;
is ($hit->name, 'Contig147');
is ($hit->description, 'Contig147.fa');
is ($hit->length, 1086);

my $hsp = $hit->next_hsp;
is ($hsp->query->start, 36);
is ($hsp->query->end, 132);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 191);
is ($hsp->hit->end, 286);
is ($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 133);
is ($hsp->query->end, 191);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 343);
is ($hsp->hit->end, 401);
is ($hsp->hit->strand, 1);

# parse align format 3
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data crypto.sim4-3))
			    );
$r = $parser->next_result;
is ($r->query_name, 'cn416');
is ($r->query_length, 630);
$hit = $r->next_hit;
is ($hit->name, 'Contig147');
is ($hit->description, 'Contig147.fa');
is ($hit->length, 1086);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 36);
is ($hsp->query->end, 132);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 191);
is ($hsp->hit->end, 286);
is ($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 133);
is ($hsp->query->end, 191);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 343);
is ($hsp->hit->end, 401);
is ($hsp->hit->strand, 1);

# parse align format 4
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data crypto.sim4-4))
			    );
$r = $parser->next_result;
is ($r->query_name, 'cn416');
is ($r->query_length, 630);

$hit = $r->next_hit;
is ($hit->name, 'Contig147');
is ($hit->length, 1086);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 36);
is ($hsp->query->end, 132);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 191);
is ($hsp->hit->end, 286);
is ($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 133);
is ($hsp->query->end, 191);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 343);
is ($hsp->hit->end, 401);
is ($hsp->hit->strand, 1);


# do the other sim4 files
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data sim4.rev))
			    );
$r = $parser->next_result;
is ($r->query_name, '/nfs/disk21/birney/prog/wise2/example/human.rev');
is ($r->query_length, 5368);
$hit = $r->next_hit;
is ($hit->name, 'HSHNCPA1');
is ($hit->description, 'temp.cdna');
is ($hit->length, 1198);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 486);
is ($hsp->query->end, 503);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 10);
is ($hsp->hit->end, 27);
is ($hsp->hit->strand, -1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 1048);
is ($hsp->query->end, 1117);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 194);
is ($hsp->hit->end, 265);
is ($hsp->hit->strand, -1);

# do the other sim4 files fwd
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data sim4.for.for))
			    );
$r = $parser->next_result;
is ($r->query_name, 'human.genomic');
is ($r->query_length, 5368);
$hit = $r->next_hit;
is ($hit->name, 'hs_est');
is ($hit->description, 'est.for');
is ($hit->length, 479);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 695);
is ($hsp->query->end, 813);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 1);
is ($hsp->hit->end, 119);
is ($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 1377);
is ($hsp->query->end, 1500);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 120);
is ($hsp->hit->end, 243);
is ($hsp->hit->strand, 1);

# do the other sim4 files fwd rev
$parser = new Bio::SearchIO(-format => 'sim4',
			    -file   => 
			    Bio::Root::IO->catfile(qw(t data sim4.for.rev))
			    );
$r = $parser->next_result;
is ($r->query_name, 'human.genomic');
is ($r->query_length, 5368);
$hit = $r->next_hit;
is ($hit->name, 'REVCOMP');
is ($hit->description, 'hn_est.rev');
is ($hit->length, 479);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 695);
is ($hsp->query->end, 813);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 1);
is ($hsp->hit->end, 119);
is ($hsp->hit->strand, -1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 1377);
is ($hsp->query->end, 1500);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 120);
is ($hsp->hit->end, 243);
is ($hsp->hit->strand, -1);
