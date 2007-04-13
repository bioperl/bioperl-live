# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
# Created: Thu Dec 14 13:57:04 GMT 2000
# By Joseph A.L. Insana, <insana@ebi.ac.uk>
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
my $error;
BEGIN {     
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    eval { require IO::String };
    if( $@ ) {
		plan skip_all => "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.";
	} else {
		plan tests => $NUMTESTS = 49;
	}
	use_ok('Bio::LiveSeq::IO::BioPerl');
	use_ok('Bio::Root::IO');
}

my $loader=Bio::LiveSeq::IO::BioPerl->load(-db=>"EMBL", 
					   -file=>Bio::Root::IO->catfile("t","data","factor7.embl"));
ok $loader;
my $gene=$loader->gene2liveseq(-gene_name => "factor7");
ok $gene;
ok ref($gene), "Bio::LiveSeq::Gene";
is $gene->name, "factor7";
is $gene->get_DNA->alphabet, "dna";
is $gene->get_DNA->display_id, "HSCFVII";
is $gene->get_DNA->accession_number, "J02933";
is $gene, $gene->get_DNA->gene;
is $gene->get_DNA->desc, "Human blood coagulation factor VII gene, complete cds.";
is $gene->get_DNA->source, "Homo sapiens";
is $gene->get_DNA->start, 1;
is $gene->get_DNA->end, 12850;
is $gene->maxtranscript->start, 487;
is $gene->maxtranscript->end, 12686;
is $gene->upbound, 487;
is $gene->downbound, 12686;
ok not(defined($gene->get_Repeat_Units));

my @exons   = @{$gene->get_Exons};
my @introns = @{$gene->get_Introns};
is scalar(@exons), 9;
is scalar(@introns), 8;
is $introns[4]->desc, "Intron D";
is $introns[4]->start, 6592;
is $introns[4]->end, 8306;
is $exons[1]->desc, "optional";
is $exons[4]->end, 6591;

my $transcript  = $gene->get_Transcripts->[0];
my $translation = $gene->get_Translations->[0];
is $transcript , $translation->get_Transcript;
is $translation , $transcript->get_Translation;

@exons = $transcript->all_Exons;
is $exons[4]->end , 6591;
is $exons[4]->length , 114;
is $transcript->upstream_seq, "tcaacaggcaggggcagcactgcagagatttcatc";
is substr($transcript->downstream_seq,0,16), "cccagcagccctggcc";
is $transcript->position($transcript->label(666)), 666;
is $transcript->position($transcript->label(666),9419), 95;
is $transcript->labelsubseq(8447,undef,9419), "gt";
is $transcript->labelsubseq(8447,2), "gt";
is $gene->get_DNA->labelsubseq(8447,2), "gg";
is substr($gene->get_DNA->labelsubseq(8447,undef,9419),0,16), "ggtgaccaggcttcat";
is $gene->get_DNA, $transcript->{seq};
my ($nothing,$whichexon) = $transcript->in_which_Exon(9419);
is $whichexon , 7;
is $transcript->frame(9419) , 1;
is $transcript->frame(9420) , 2;
is substr($translation->seq,0,16), "MVSQALRLLCLLLGLQ";
is substr($transcript->seq,0,32), "atggtctcccaggccctcaggctcctctgcct";
ok $transcript->translation_table(2);
is $transcript->translation_table , 2;
is substr($translation->seq,0,16), "MVSQAL*"; # mitochondrial table creates stop codon
is $gene->verbose(2), 2;
ok $gene->delete_Obj(); # to free all memory, deleting circular references

