#!/usr/local/bin/perl
#
# Script for checking and demonstrating the Genscan parser module.
#
# usage: pipe the Genscan prediction result to this script. Optionally,
# specify the sequence file (FASTA format) fed to Genscan by --seq myseq.tfa.
#
# BTW you can use this script for testing the MZEF module as well with a few
# modifications: replace every occurrence of Bio::Tools::Genscan with
# Bio::Tools::MZEF.
#
# Hilmar Lapp <hlapp@gmx.net>
#

use Getopt::Long;
use Bio::Tools::Genscan;
use Bio::Tools::Prediction::Gene;
use Bio::SeqIO;

my $in = \*STDIN;
my $out = \*STDOUT;

$opt_ok = GetOptions("seq=s", \$seqfile);

my $seq;

if(defined($seqfile)) {
    my $seqin = Bio::SeqIO->new('-file' => $seqfile, '-format' => "fasta");
    $seq = $seqin->next_seq();
    $seqin->close();
}
my $seqout = Bio::SeqIO->new('-fh' => $out);

my $genscan = Bio::Tools::Genscan->new('-fh' => $in);

while($gene = $genscan->next_prediction()) {
    $gene->attach_seq($seq) if $seq;
    
    my $fea;
    
    print $out "\n======= Next prediction ==========\n\n";
    foreach $fea ($gene->exons()) {
	print $out "EXON ", $fea->gff_string(), "\n";
    }
    foreach $fea ($gene->promotors()) {
	print $out "PROM ", $fea->gff_string(), "\n";
    }
    foreach $fea ($gene->poly_A_sites()) {
	print $out "POLA ", $fea->gff_string(), "\n";
    }
    # sequences
    my $cds;
    my $outseq;

    if($gene->predicted_protein()) {
	$seqout->write_seq($gene->predicted_protein());
    }
    if($seq) {
	$cds = Bio::PrimarySeq->new('-seq' => $gene->cds(1),
				    '-id' => "cds_" . $gene->primary_tag());
	$outseq = $cds->translate();
	$seqout->write_seq($outseq);
    }
    if($gene->predicted_cds()) {
	$seqout->write_seq($gene->predicted_cds());
    }
    if($seq) {
	$seqout->write_seq($cds);
    }	
}

