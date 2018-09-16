#!/usr/bin/perl

=head1 NAME

reverse-translate.pl

=head1 DESCRIPTION

Reverse-translates a nucleotide sequence using the most frequent codons.
Requires an input sequence file and a nucleotide sequence file containing 
one sequence comprised of one or more ORFs. This file supplies the codon
frequency data and will be parsed starting at the first triplet in the sequence.

=head1 OPTIONS

  -i   Input sequence, amino acid
  -c   Input sequence, nucleotide ORFs

Example:

 reverse-translate.pl -i ~/bioperl-live/t/data/cysprot.fa -c ~/bioperl-live/t/data/HUMBETGLOA.fa 

=cut

use strict;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Bio::Tools::SeqStats;
use Bio::CodonUsage::Table;
use Getopt::Long;

my ($codonFile,$seqFile);

GetOptions( "c=s" => \$codonFile,
            "i=s" => \$seqFile );

die "Need input sequence and file containing coding regions"
  if ( !$codonFile || !$seqFile );

my $codonIn = Bio::SeqIO->new(-file => $codonFile,
										-format => 'fasta');

my $codonSeq = $codonIn->next_seq;

my $codonStats = Bio::Tools::SeqStats->count_codons($codonSeq);

my $codonUsage = Bio::CodonUsage::Table->new(-data => $codonStats );

my $codonTable = Bio::Tools::CodonTable->new;

my $seqIn = Bio::SeqIO->new(-file => $seqFile);

my $seq = $seqIn->next_seq;

my $rvSeq = $codonTable->reverse_translate_best($seq,$codonUsage);

print $rvSeq,"\n";

__END__
