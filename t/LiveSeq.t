# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
# Created: Thu Dec 14 13:57:04 GMT 2000
# By Joseph A.L. Insana, <insana@ebi.ac.uk>
#
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
BEGIN { $| = 1; print "1..44\n";
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::LiveSeq::IO::BioPerl;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

my $loader=Bio::LiveSeq::IO::BioPerl->load(-db=>"EMBL", -file=>"factor7.embl");
test 2, $loader;
my $gene=$loader->gene2liveseq(-gene_name => "factor7");
test 3, $gene->name eq "factor7";
test 4, $gene->get_DNA->moltype eq "dna";
test 5, $gene->get_DNA->display_id eq "HSCFVII";
test 6, $gene->get_DNA->accession_number eq "J02933";
test 7, $gene == $gene->get_DNA->gene;
test 8, $gene->get_DNA->description eq "Human blood coagulation factor VII gene, complete cds.";
test 9, $gene->get_DNA->source eq "Homo sapiens";
test 10, $gene->get_DNA->start == 1;
test 11, $gene->get_DNA->end == 12850;
test 12, $gene->maxtranscript->start == 487;
test 13, $gene->maxtranscript->end == 12686;
test 14, $gene->upbound == 487;
test 15, $gene->downbound == 12686;
test 16, not(defined($gene->get_Repeat_Units));
my @exons=@{$gene->get_Exons};
my @introns=@{$gene->get_Introns};
test 17, scalar(@exons) == 9;
test 18, scalar(@introns) == 8;
test 19, $introns[4]->description eq "Intron D";
test 20, $introns[4]->start == 6592;
test 21, $introns[4]->end == 8306;
test 22, $exons[1]->description eq "optional";
test 23, $exons[4]->end == 6591;
my $transcript=$gene->get_Transcripts->[0];
my $translation=$gene->get_Translations->[0];
test 24, $transcript == $translation->get_Transcript;
test 25, $translation == $transcript->get_Translation;
@exons=$transcript->all_Exons;
test 26, $exons[4]->end == 6591;
test 27, $exons[4]->length == 114;
test 28, $transcript->upstream_seq eq "tcaacaggcaggggcagcactgcagagatttcatc";
test 29, substr($transcript->downstream_seq,0,16) eq "cccagcagccctggcc";
test 30, $transcript->position($transcript->label(666)) eq 666;
test 31, $transcript->position($transcript->label(666),9419) eq 95;
test 32, $transcript->labelsubseq(8447,undef,9419) eq "gt";
test 33, $transcript->labelsubseq(8447,2) eq "gt";
test 34, $gene->get_DNA->labelsubseq(8447,2) eq "gg";
test 35, substr($gene->get_DNA->labelsubseq(8447,undef,9419),0,16) eq "ggtgaccaggcttcat";
test 36, $gene->get_DNA eq $transcript->{seq};
my (undef,$whichexon)=$transcript->in_which_Exon(9419);
test 37, $whichexon == 7;
test 39, $transcript->frame(9419) == 1;
test 40, $transcript->frame(9420) == 2;
test 41, substr($translation->seq,0,16) eq "MVSQALRLLCLLLGLQ";
test 42, substr($transcript->seq,0,32) eq "atggtctcccaggccctcaggctcctctgcct";
test 43, $transcript->translation_table(2) && $transcript->translation_table == 2;
test 44, substr($translation->seq,0,16) eq "MVSQAL*"; # mitochondrial table creates stop codon
