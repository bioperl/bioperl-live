# -*-Perl-*- Test Harness script for Bioperl
# $Id: Genewise.t 11733 2007-10-26 18:22:10Z jason $

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 33);
	
    use_ok('Bio::Tools::Genewise');
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
my ($phase) = $e[0]->get_tag_values('phase');
is ($phase,0);
my ($sf)= $e[0]->get_tag_values('supporting_feature');
is ($sf->feature1->seq_id,'Scaffold_2042.1');
is ($sf->feature1->start,22265);
is ($sf->feature1->end,22396);
is ($sf->feature2->seq_id,'SINFRUP00000067802');
is ($sf->feature2->start,1);
is ($sf->feature2->end,44);
is ($sf->feature1->end,22396);

open my $FH, '<', $inputfilename or die "Could not read file '$inputfilename': $!\n";
$parser = Bio::Tools::Genewise->new(-fh => $FH);
while (my $gene = $parser->next_prediction){
    push @gene, $gene;
}
@t = $gene[0]->transcripts;
@e = $t[0]->exons;

is (scalar($t[0]->exons), 18);
is ($t[0]->start, 22265);
is ($t[0]->end, 37062);
is ($e[0]->start,22265);
is ($e[0]->end, 22396);
($phase) = $e[0]->get_tag_values('phase');
is ($phase,0);
($sf)= $e[0]->get_tag_values('supporting_feature');
is ($sf->feature1->seq_id,'Scaffold_2042.1');
is ($sf->feature1->start,22265);
is ($sf->feature1->end,22396);
is ($sf->feature2->seq_id,'SINFRUP00000067802');
is ($sf->feature2->start,1);
is ($sf->feature2->end,44);
is ($sf->feature1->end,22396);
