# -*-Perl-*-
## Bioperl Test Harness Script for Modules

use strict;

BEGIN
{
    eval { require Test; };
    if ($@)
    {
        use lib 't';
    }

    use Test;
    use vars qw($NTESTS);
    $NTESTS = 26;
    plan tests => $NTESTS;
}

use Bio::Tools::Geneid;
use Bio::SeqIO;
ok(1);

my $inputfilename = Bio::Root::IO->catfile("t", "data", "geneid_1.0.out");
my $parser = Bio::Tools::Geneid->new(-file => $inputfilename);
my @genes;

while (my $gene= $parser->next_prediction)
{
    push(@genes, $gene);
}

my @transcripts = $genes[0]->transcripts;
my @exons = $transcripts[0]->exons;

ok($transcripts[0]->seq_id, '10');
ok($exons[0]->seq_id, '10');
ok($transcripts[0]->source_tag, 'geneid');
ok($exons[0]->source_tag, 'geneid');
ok($transcripts[0]->primary_tag, 'transcript');
ok($exons[0]->primary_tag, 'Initial');

ok(scalar($transcripts[0]->exons), 2);
ok($transcripts[0]->start, 6090);
ok($transcripts[0]->end, 7276);
ok($transcripts[0]->score, 36.87);
ok($transcripts[0]->strand, 1);
ok($exons[0]->start, 6090);
ok($exons[0]->end, 6155);
ok($exons[0]->score, '1.40');
ok($exons[0]->strand, 1);

my ($type) = $exons[0]->get_tag_values('Type');
ok($type, 'Initial');

my ($phase) = $exons[0]->get_tag_values('phase');
ok($phase, 0);

my ($end_phase) = $exons[0]->get_tag_values('end_phase');
ok($end_phase, 0);

my ($start_signal_score) = $exons[0]->get_tag_values('start_signal_score');
ok($start_signal_score, 2.15);

my ($end_signal_score) = $exons[0]->get_tag_values('end_signal_score');
ok($end_signal_score, 3.63);

my ($coding_potential_score) = $exons[0]->get_tag_values('coding_potential_score');
ok($coding_potential_score, 12.34);

my ($homology_score) = $exons[0]->get_tag_values('homology_score');
ok($homology_score, '0.00');

ok(scalar(@genes), 3);

@transcripts = $genes[1]->transcripts;
ok(scalar($transcripts[0]->exons), 5);

@transcripts = $genes[2]->transcripts;
ok(scalar($transcripts[0]->exons), 1);
