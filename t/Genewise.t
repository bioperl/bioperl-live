# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 53);
	
    use_ok('Bio::Tools::Genewise');
    use_ok('Bio::SearchIO');
}

my $inputfilename= test_input_file('genewise.out');
my $parser = Bio::Tools::Genewise->new(-file => $inputfilename);
my @gene;
while (my $gene= $parser->next_prediction){
    push @gene, $gene;
}
my @t = $gene[0]->transcripts;
my @e = $t[0]->exons;

is ($t[0]->seq_id, 'Scaffold_2042.1');
is ($e[0]->seq_id, 'Scaffold_2042.1');
is ($t[0]->source_tag, 'genewise');
is ($e[0]->source_tag, 'genewise');
is ($t[0]->primary_tag, 'transcript');
is ($e[0]->primary_tag, 'exon');

is (scalar($t[0]->exons), 18);
is ($t[0]->start, 22265);
is ($t[0]->end, 37062);
is ($e[0]->start,22265);
is ($e[0]->end, 22396);
my ($phase) = $e[0]->each_tag_value('phase');
is ($phase,0);
my ($sf)= $e[0]->each_tag_value('supporting_feature');
is ($sf->feature1->seq_id,'Scaffold_2042.1');
is ($sf->feature1->start,22265);
is ($sf->feature1->end,22396);
is ($sf->feature2->seq_id,'SINFRUP00000067802');
is ($sf->feature2->start,1);
is ($sf->feature2->end,44);
is ($sf->feature1->end,22396);

open(FH,$inputfilename);
$parser = Bio::Tools::Genewise->new(-fh=>\*FH);
while (my $gene= $parser->next_prediction){
    push @gene, $gene;
}
@t = $gene[0]->transcripts;
@e = $t[0]->exons;

is (scalar($t[0]->exons), 18);
is ($t[0]->start, 22265);
is ($t[0]->end, 37062);
is ($e[0]->start,22265);
is ($e[0]->end, 22396);
($phase) = $e[0]->each_tag_value('phase');
is ($phase,0);
($sf)= $e[0]->each_tag_value('supporting_feature');
is ($sf->feature1->seq_id,'Scaffold_2042.1');
is ($sf->feature1->start,22265);
is ($sf->feature1->end,22396);
is ($sf->feature2->seq_id,'SINFRUP00000067802');
is ($sf->feature2->start,1);
is ($sf->feature2->end,44);
is ($sf->feature1->end,22396);

$parser = Bio::SearchIO->new(-file => test_input_file('genewise.out'),
			    -format   => 'wise',
			    -wisetype => 'genewise');
my $result = $parser->next_result;
my $hit = $result->next_hit;
is($result->query_name, 'SINFRUP00000067802');
is($hit->name, 'Scaffold_2042.1');

is($hit->score, 2054.68);
my $hsp = $hit->next_hsp;

is($hsp->query->start,22265);
is($hsp->query->end,22396);
is($hsp->query->strand,1);
is($hsp->query->score, 2054.68);

is($hsp->hit->start,1);
is($hsp->hit->end,44);
is($hsp->hit->strand,0);
is($hsp->hit->score, 2054.68);

$hsp = $hit->next_hsp;

is($hsp->query->start,24224);
is($hsp->query->end,24328);

is($hsp->hit->start,45);
is($hsp->hit->end,79);

$hsp = $hit->next_hsp;

is($hsp->query->start,24471);
is($hsp->query->end,24513);

is($hsp->hit->start,80);
is($hsp->hit->end,93);
