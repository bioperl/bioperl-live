# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

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
    plan tests => 16; 
}

if( $error == 1 ) {
    exit(0);
}

my $debug = -1;

use Bio::Align::Utilities qw(aa_to_dna_aln bootstrap_replicates);
use Bio::AlignIO;
use Bio::Root::IO;
use Bio::SeqIO;

my $in = new Bio::AlignIO(-format => 'clustalw',
			  -file   => Bio::Root::IO->catfile
			  ('t','data','pep-266.aln'));
my $aln = $in->next_aln();
ok($aln);
$in->close();

my $seqin = new Bio::SeqIO(-format => 'fasta',
			   -file   => Bio::Root::IO->catfile
			   ('t','data','cds-266.fas'));
# get the cds sequences
my %cds_seq;
while( my $seq = $seqin->next_seq ) {
    $cds_seq{$seq->display_id} = $seq;
}

my $cds_aln = &aa_to_dna_aln($aln,\%cds_seq);

my @aa_seqs = $aln->each_seq;

for my $cdsseq ( $cds_aln->each_seq ) {
    my $peptrans = $cdsseq->translate();
    my $aaseq = shift @aa_seqs;
    ok($peptrans->seq(),$aaseq->seq());
}

my $bootstraps = &bootstrap_replicates($aln,10);

ok(scalar @$bootstraps, 10);
