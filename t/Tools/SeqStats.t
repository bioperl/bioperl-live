# -*-Perl-*- Test Harness script for Bioperl
# $Id: SeqStats.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN { 
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 7);
	
	use_ok('Bio::SeqIO');
	use_ok('Bio::Tools::SeqStats');
}

my $str = Bio::SeqIO->new(-file=> test_input_file('multifa.seq'), '-format' => 'Fasta');
my $seqobj = $str->next_seq();
ok defined $seqobj, 'new Bio::Root::IO object';

my $seq_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);
ok defined $seq_stats && $seq_stats, 'new Bio::Tools:SeqStats object';

# eg for DNA sequence
my $hash_ref = $seq_stats->count_monomers();  
is ( $hash_ref->{'A'}, 80 , 'count_monomers()');

$hash_ref = $seq_stats-> count_codons();  
ok defined $hash_ref && $hash_ref , 'count_codons()';

my $weight = $seq_stats->get_mol_wt();
ok defined $weight && $weight , 'get_mol_wt()' ;
