# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 8; 
}

use Bio::SeqIO;
use Bio::Tools::SeqWords;
use Bio::Tools::SeqStats;
use Bio::Root::IO;

ok(1);
my $str = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","multifa.seq"), '-format' => 'Fasta');
my $seqobj= $str->next_seq();
ok $seqobj;

my $words = Bio::Tools::SeqWords->new('-seq' => $seqobj);
my $hash = $words->count_words(6);
ok ($words);
ok ($hash);

my $seq_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);

ok $seq_stats;

my $hash_ref = $seq_stats->count_monomers();  # eg for DNA sequence

ok ( $hash_ref->{'A'}, 80 );

$hash_ref = $seq_stats-> count_codons();  

ok $hash_ref;

my $weight = $seq_stats->get_mol_wt();
ok $weight;

