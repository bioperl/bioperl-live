# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    plan tests => 11 ; 
}

use_ok('Bio::SeqIO');
use_ok('Bio::Tools::SeqWords');
use_ok('Bio::Tools::SeqStats');
use_ok('Bio::Root::IO');

my $str = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","multifa.seq"), '-format' => 'Fasta');
my $seqobj= $str->next_seq();
ok defined $seqobj, 'new Bio::Root::IO object';

my $words = Bio::Tools::SeqWords->new('-seq' => $seqobj);
my $hash = $words->count_words(6);
ok (defined $words, 'new Bio::Tools::SeqWords object');
ok (defined $hash, 'count_words');

my $seq_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);
ok defined $seq_stats && $seq_stats, 'new Bio::Tools:SeqStats object';

# eg for DNA sequence
my $hash_ref = $seq_stats->count_monomers();  
is ( $hash_ref->{'A'}, 80 , 'count_monomers()');

$hash_ref = $seq_stats-> count_codons();  
ok defined $hash_ref && $hash_ref , 'count_codons()';

my $weight = $seq_stats->get_mol_wt();
ok defined $weight && $weight , 'get_mol_wt()' ;

