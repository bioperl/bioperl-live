# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {
    eval { require Test::More; };
    if( $@ ) {
		use lib 't/lib';
    }

    use Test::More;
    plan tests => 20;
	use_ok('Bio::Align::Utilities',qw(aa_to_dna_aln bootstrap_replicates) );
	use_ok('Bio::AlignIO');
	use_ok('Bio::Root::IO');
	use_ok('Bio::SeqIO');
}

my $debug = -1;

my $in = Bio::AlignIO->new(-format => 'clustalw',
			  -file   => Bio::Root::IO->catfile
			  ('t','data','pep-266.aln'));
my $aln = $in->next_aln();
isa_ok($aln, 'Bio::Align::AlignI');
$in->close();

my $seqin = Bio::SeqIO->new(-format => 'fasta',
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
    is($peptrans->seq(),$aaseq->seq());
}

my $bootstraps = &bootstrap_replicates($aln,10);

is(scalar @$bootstraps, 10);
