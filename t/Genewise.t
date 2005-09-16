# -*-Perl-*-
## Bioperl Test Harness Script for Modules


use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 51;
    plan tests => $NTESTS;
}
use Bio::Tools::Genewise;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Root::IO;

END {
	for ( $Test::ntest..$NTESTS ) {
		skip("Cannot run remaining Genewise tests, skipping.",1);
	}
}

my $inputfilename= Bio::Root::IO->catfile("t","data","genewise.out");
my $parser = Bio::Tools::Genewise->new(-file => $inputfilename);
my @gene;
while (my $gene= $parser->next_prediction){
    push @gene, $gene;
}
my @t = $gene[0]->transcripts;
my @e = $t[0]->exons;

ok ($t[0]->seq_id, 'Scaffold_2042.1');
ok ($e[0]->seq_id, 'Scaffold_2042.1');
ok ($t[0]->source_tag, 'genewise');
ok ($e[0]->source_tag, 'genewise');
ok ($t[0]->primary_tag, 'transcript');
ok ($e[0]->primary_tag, 'exon');

ok (scalar($t[0]->exons), 18);
ok ($t[0]->start, 22265);
ok ($t[0]->end, 37062);
ok ($e[0]->start,22265);
ok ($e[0]->end, 22396);
my ($phase) = $e[0]->each_tag_value('phase');
ok ($phase,0);
my ($sf)= $e[0]->each_tag_value('supporting_feature');
ok ($sf->feature1->seq_id,'Scaffold_2042.1');
ok ($sf->feature1->start,22265);
ok ($sf->feature1->end,22396);
ok ($sf->feature2->seq_id,'SINFRUP00000067802');
ok ($sf->feature2->start,1);
ok ($sf->feature2->end,44);
ok ($sf->feature1->end,22396);

open(FH,$inputfilename);
$parser = Bio::Tools::Genewise->new(-fh=>\*FH);
while (my $gene= $parser->next_prediction){
    push @gene, $gene;
}
@t = $gene[0]->transcripts;
@e = $t[0]->exons;

ok (scalar($t[0]->exons), 18);
ok ($t[0]->start, 22265);
ok ($t[0]->end, 37062);
ok ($e[0]->start,22265);
ok ($e[0]->end, 22396);
($phase) = $e[0]->each_tag_value('phase');
ok ($phase,0);
($sf)= $e[0]->each_tag_value('supporting_feature');
ok ($sf->feature1->seq_id,'Scaffold_2042.1');
ok ($sf->feature1->start,22265);
ok ($sf->feature1->end,22396);
ok ($sf->feature2->seq_id,'SINFRUP00000067802');
ok ($sf->feature2->start,1);
ok ($sf->feature2->end,44);
ok ($sf->feature1->end,22396);

$parser = new Bio::SearchIO(-file => 
			    Bio::Root::IO->catfile(qw(t data genewise.out)),
			    -format   => 'wise',
			    -wisetype => 'genewise');
my $result = $parser->next_result;
skip(1,'swapping query/name need to reconsider how this done');
#ok($result->query_name, 'SINFRUP00000067802');
my $hit = $result->next_hit;
skip(1,'swapping query/name need to reconsider how this done');
#ok($hit->name, 'Scaffold_2042.1');
ok($hit->score, 2054.68);
my $hsp = $hit->next_hsp;

ok($hsp->query->start,22265);
ok($hsp->query->end,22396);
ok($hsp->query->strand,1);
ok($hsp->query->score, 2054.68);

ok($hsp->hit->start,1);
ok($hsp->hit->end,44);
ok($hsp->hit->strand,0);
ok($hsp->hit->score, 2054.68);

$hsp = $hit->next_hsp;

ok($hsp->query->start,24224);
ok($hsp->query->end,24328);

ok($hsp->hit->start,45);
ok($hsp->hit->end,79);

$hsp = $hit->next_hsp;

ok($hsp->query->start,24471);
ok($hsp->query->end,24513);

ok($hsp->hit->start,80);
ok($hsp->hit->end,93);
