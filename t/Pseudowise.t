#!/usr/local/bin/perl
# -*-Perl-*-
## Bioperl Test Harness Script for Modules

use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 12;
    plan tests => $NTESTS;
}
use Bio::Tools::Run::Pseudowise;
use Bio::Root::IO;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("pseudowise program not found. Skipping. (Be sure you have the wise package > 2.2.0)",1);
    }
}

ok(1);
my $verbose = -1;
my @params = ('dymem', 'linear', 'kbyte', '5000');
my  $factory = Bio::Tools::Run::Pseudowise->new(@params);
ok $factory->isa('Bio::Tools::Run::Pseudowise');
my $bequiet = 1;
$factory->quiet($bequiet);  # Suppress pseudowise messages to terminal


#test with one file with 2 sequences
my $inputfilename = Bio::Root::IO->catfile("t","data","ps1.fa");
my $seq1 = Bio::Seq->new();
my $seq2 = Bio::Seq->new();
my $seq3 = Bio::Seq->new();
my $seqstream = Bio::SeqIO->new(-file => $inputfilename, -fmt => 'Fasta');
$seq1 = $seqstream->next_seq();
$seq2 = $seqstream->next_seq();
$seq3 = $seqstream->next_seq();


my $pseudowise_present = $factory->exists_pseudowise();
unless ($pseudowise_present) {
    warn("Pseudowise program not found. Skipping tests $Test::ntest to $NTESTS.\n");
    exit 0;
}
my @feat = $factory->predict_genes($seq1, $seq2, $seq3);
my $geneno = scalar(@feat);
my @subfeat = $feat[0]->sub_SeqFeature();
my $exonno = scalar(@subfeat);

ok($geneno, 1);
ok($exonno, 1);
ok($feat[0]->isa("Bio::SeqFeatureI"));
ok($subfeat[0]->isa("Bio::SeqFeatureI"));
ok($feat[0]->primary_tag, 'pseudogene');
ok($subfeat[0]->primary_tag, 'exon');
ok($feat[0]->start, 865);
ok($subfeat[0]->start, 865);
ok($feat[0]->end, 897);
ok($subfeat[0]->end, 897);






