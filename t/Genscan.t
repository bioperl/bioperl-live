#!/usr/local/bin/perl
# -*-Perl-*-
# ## Bioperl Test Harness Script for Modules
#
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
use Bio::Tools::Run::Genscan;
use Bio::Root::IO;

ok(1);
my $verbose = -1;

my $paramfile = ("data/HumanIso.smat");
my  $factory = Bio::Tools::Run::Genscan->new(MATRIX=>$paramfile);
ok $factory->isa('Bio::Tools::Run::Genscan');


#test with one file with 2 sequences
my $inputfilename = Bio::Root::IO->catfile("data","Genscan.FastA");
#my $inputfilename = Bio::Root::IO->catfile("code1.FastA");
my $seq1 = Bio::Seq->new();
my $seqstream = Bio::SeqIO->new(-file => $inputfilename, -fmt => 'Fasta');
$seq1 = $seqstream->next_seq();


my $genscan_present = $factory->exists_genscan();
unless ($genscan_present) {
        warn("Genscan program not found. Skipping tests $Test::ntest to $NTESTS.\n");
            exit 0;
}
my @feat = $factory->predict_genes($seq1); 
my $protein = $feat[0]->predicted_protein();
ok($feat[0]->isa("Bio::SeqFeatureI"));
ok($protein->isa("Bio::PrimarySeqI"));
