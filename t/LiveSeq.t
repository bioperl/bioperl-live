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
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $NUMTESTS = 48;
    plan tests => $NUMTESTS;
    $error = 0;
    eval { require 'IO/String.pm' };
    if( $@ ) {
	print STDERR "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n";
	for( 1..$NUMTESTS ) {
	    skip(1,"IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests");
	}
	$error = 1; 
    }
}

if( $error ==  1 ) {
    exit(0);
}

require Bio::LiveSeq::IO::BioPerl;
require Bio::Root::IO;

ok(1);

my $loader=Bio::LiveSeq::IO::BioPerl->load(-db=>"EMBL", 
					   -file=>Bio::Root::IO->catfile("t","data","factor7.embl"));
ok $loader;
my $gene=$loader->gene2liveseq(-gene_name => "factor7");
ok $gene;
ok ref($gene), "Bio::LiveSeq::Gene";
ok $gene->name, "factor7";
ok $gene->get_DNA->moltype, "dna";
ok $gene->get_DNA->display_id, "HSCFVII";
ok $gene->get_DNA->accession_number, "J02933";
ok $gene, $gene->get_DNA->gene;
ok $gene->get_DNA->desc, "Human blood coagulation factor VII gene, complete cds.";
ok $gene->get_DNA->source, "Homo sapiens";
ok $gene->get_DNA->start, 1;
ok $gene->get_DNA->end, 12850;
ok $gene->maxtranscript->start, 487;
ok $gene->maxtranscript->end, 12686;
ok $gene->upbound, 487;
ok $gene->downbound, 12686;
ok not(defined($gene->get_Repeat_Units));

my @exons   = @{$gene->get_Exons};
my @introns = @{$gene->get_Introns};
ok scalar(@exons), 9;
ok scalar(@introns), 8;
ok $introns[4]->desc, "Intron D";
ok $introns[4]->start, 6592;
ok $introns[4]->end, 8306;
ok $exons[1]->desc, "optional";
ok $exons[4]->end, 6591;

my $transcript  = $gene->get_Transcripts->[0];
my $translation = $gene->get_Translations->[0];
ok $transcript , $translation->get_Transcript;
ok $translation , $transcript->get_Translation;

@exons = $transcript->all_Exons;
ok $exons[4]->end , 6591;
ok $exons[4]->length , 114;
ok $transcript->upstream_seq, "tcaacaggcaggggcagcactgcagagatttcatc";
ok substr($transcript->downstream_seq,0,16), "cccagcagccctggcc";
ok $transcript->position($transcript->label(666)), 666;
ok $transcript->position($transcript->label(666),9419), 95;
ok $transcript->labelsubseq(8447,undef,9419), "gt";
ok $transcript->labelsubseq(8447,2), "gt";
ok $gene->get_DNA->labelsubseq(8447,2), "gg";
ok substr($gene->get_DNA->labelsubseq(8447,undef,9419),0,16), "ggtgaccaggcttcat";
ok $gene->get_DNA, $transcript->{seq};
my ($nothing,$whichexon) = $transcript->in_which_Exon(9419);
ok $whichexon , 7;
ok $transcript->frame(9419) , 1;
ok $transcript->frame(9420) , 2;
ok substr($translation->seq,0,16), "MVSQALRLLCLLLGLQ";
ok substr($transcript->seq,0,32), "atggtctcccaggccctcaggctcctctgcct";
ok $transcript->translation_table(2);
ok $transcript->translation_table , 2;
ok substr($translation->seq,0,16), "MVSQAL*"; # mitochondrial table creates stop codon
ok $gene->verbose(2), 2;
ok $gene->delete_Obj(); # to free all memory, deleting circular references

