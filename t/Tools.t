# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..6\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::SeqIO;
use Bio::Tools::SeqWords;
use Bio::Tools::SeqStats;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}


$str = Bio::SeqIO->new(-file=> 't/multifa.seq', '-format' => 'Fasta');
$seqobj= $str->next_seq();

$words = Bio::Tools::SeqWords->new($seqobj);
$hash = $words->count_words(6);
$hash = undef; # shut warn up
test 2,1;

$seq_stats  =  Bio::Tools::SeqStats->new($seqobj);

test 3,1;

$hash_ref = $seq_stats->count_monomers();  # eg for DNA sequence

test 4, ( $hash_ref->{'A'} == 80 );

$hash_ref = $seq_stats-> count_codons();  

test 5, 1;

$weight = $seq_stats->get_mol_wt();
$weight = 0; # shut warn up
test 6, 1;
