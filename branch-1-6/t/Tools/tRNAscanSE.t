# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14);
	
	use_ok('Bio::Tools::tRNAscanSE');
}

my $verbose = test_debug();

my $parser = Bio::Tools::tRNAscanSE->new(-verbose => $verbose,
					 -file => test_input_file('yeast.tRNAscanSE'));

isa_ok($parser, 'Bio::Tools::tRNAscanSE') ;

my @genes;
while( my $gene = $parser->next_prediction ) {
    push @genes, $gene;
}

is (scalar(@genes), 287);
is($genes[2]->seq_id, 'I', 'seq_id');
my ($codon) = $genes[2]->get_tag_values('Codon');
is($codon, 'TTG', 'codon');
is($genes[2]->start, 181135, 'start');
is($genes[2]->end, 181248, 'end');
is($genes[2]->strand, 1, 'strand');

my @exons = $genes[2]->get_SeqFeatures ;
is ( scalar(@exons), 2, 'exons' );
is($exons[0]->end,181172, 'end' ); 
is($exons[0]->start,$genes[2]->start, 'start'); 
is($exons[1]->start,181205, 'start'); 
is($exons[1]->end,$genes[2]->end, 'end'); 
is($exons[0]->seq_id, $genes[2]->seq_id, 'seq_id');
