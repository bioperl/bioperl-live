# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;

BEGIN { plan tests => 8 }

use Bio::SeqIO;
use Bio::Tools::SeqWords;
use Bio::Tools::SeqStats;

ok(1);
my $str = Bio::SeqIO->new(-file=> 't/multifa.seq', '-format' => 'Fasta');
my $seqobj= $str->next_seq();
ok $seqobj;

my $words = Bio::Tools::SeqWords->new($seqobj);
my $hash = $words->count_words(6);
ok ($words);
ok ($hash);

my $seq_stats  =  Bio::Tools::SeqStats->new($seqobj);

ok $seq_stats;

my $hash_ref = $seq_stats->count_monomers();  # eg for DNA sequence

ok ( $hash_ref->{'A'}, 80 );

$hash_ref = $seq_stats-> count_codons();  

ok $hash_ref;

my $weight = $seq_stats->get_mol_wt();
ok $weight;
