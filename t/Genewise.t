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
    $NTESTS = 32;
    plan tests => $NTESTS;
}
use Bio::Tools::Genewise;
use Bio::SeqIO;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("Genewise program not found. Skipping. (Be sure you have the wise package )",1);
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










