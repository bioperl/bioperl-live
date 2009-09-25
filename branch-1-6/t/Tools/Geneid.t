# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 26);
	
    use_ok('Bio::Tools::Geneid');
}

my $inputfilename = test_input_file('geneid_1.0.out');
my $parser = Bio::Tools::Geneid->new(-file => $inputfilename);
my @genes;

while (my $gene= $parser->next_prediction)
{
    push(@genes, $gene);
}

my @transcripts = $genes[0]->transcripts;
my @exons = $transcripts[0]->exons;

is($transcripts[0]->seq_id, '10');
is($exons[0]->seq_id, '10');
is($transcripts[0]->source_tag, 'geneid');
is($exons[0]->source_tag, 'geneid');
is($transcripts[0]->primary_tag, 'transcript');
is($exons[0]->primary_tag, 'Initial');

is(scalar($transcripts[0]->exons), 2);
is($transcripts[0]->start, 6090);
is($transcripts[0]->end, 7276);
is($transcripts[0]->score, 36.87);
is($transcripts[0]->strand, 1);
is($exons[0]->start, 6090);
is($exons[0]->end, 6155);
is($exons[0]->score, '1.40');
is($exons[0]->strand, 1);

my ($type) = $exons[0]->get_tag_values('Type');
is($type, 'Initial');

my ($phase) = $exons[0]->get_tag_values('phase');
is($phase, 0);

my ($end_phase) = $exons[0]->get_tag_values('end_phase');
is($end_phase, 0);

my ($start_signal_score) = $exons[0]->get_tag_values('start_signal_score');
is($start_signal_score, 2.15);

my ($end_signal_score) = $exons[0]->get_tag_values('end_signal_score');
is($end_signal_score, 3.63);

my ($coding_potential_score) = $exons[0]->get_tag_values('coding_potential_score');
is($coding_potential_score, 12.34);

my ($homology_score) = $exons[0]->get_tag_values('homology_score');
is($homology_score, '0.00');

is(scalar(@genes), 3);

@transcripts = $genes[1]->transcripts;
is(scalar($transcripts[0]->exons), 5);

@transcripts = $genes[2]->transcripts;
is(scalar($transcripts[0]->exons), 1);
