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
    $NTESTS = 20;
    plan tests => $NTESTS;
}
use Bio::Tools::Genomewise;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("Genewise program not found. Skipping. (Be sure you have the wise package )",1);
    }
}


my $inputfilename= Bio::Root::IO->catfile("t","data","genomewise.out");
my $parser = Bio::Tools::Genomewise->new(-file => $inputfilename);
my @gene;
while (my $gene= $parser->next_prediction){
    push @gene, $gene;
}
my @t = $gene[0]->transcripts;
my @e = $t[0]->exons;

ok ($t[0]->source_tag, 'genomewise');
ok ($e[0]->source_tag, 'genomewise');
ok ($t[0]->primary_tag, 'transcript');
ok ($e[0]->primary_tag, 'exon');

ok (scalar($t[0]->exons), 5);
ok ($t[0]->start, 4761);
ok ($t[0]->end, 6713);
ok ($e[0]->start,4761);
ok ($e[0]->end, 4874);
my ($phase) = $e[0]->each_tag_value('phase');
ok ($phase,0);

open(FH,$inputfilename);
$parser = Bio::Tools::Genomewise->new(-fh=>\*FH);
while (my $gene= $parser->next_prediction){
    push @gene, $gene;
}
@t = $gene[1]->transcripts;
@e = $t[0]->exons;

ok ($t[0]->source_tag, 'genomewise');
ok ($e[0]->source_tag, 'genomewise');
ok ($t[0]->primary_tag, 'transcript');
ok ($e[0]->primary_tag, 'exon');

ok (scalar($t[0]->exons), 3);
ok ($t[0]->start, 9862);
ok ($t[0]->end, 10316);
ok ($e[1]->start,10024);
ok ($e[1]->end, 10211);

($phase) = $e[2]->each_tag_value('phase');
ok ($phase,2);










