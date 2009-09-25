# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 33);
	
	use_ok('Bio::Align::Utilities',qw(aa_to_dna_aln bootstrap_replicates cat) );
	use_ok('Bio::AlignIO');
	use_ok('Bio::SeqIO');
}

my $in = Bio::AlignIO->new(-format => 'clustalw',
			  -file   => test_input_file('pep-266.aln'));
my $aln = $in->next_aln();
isa_ok($aln, 'Bio::Align::AlignI');
$in->close();

my $seqin = Bio::SeqIO->new(-format => 'fasta',
			   -file   => test_input_file('cds-266.fas'));
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

my $sub_aln1=$aln->slice(1,100);
my $sub_aln2=$aln->slice(101,200);
my $sub_aln3=$aln->slice(1,200);
my $cat_aln=cat($sub_aln1, $sub_aln2);
my @seq=$sub_aln3->each_seq;
for my $seq ($cat_aln->each_seq) {
    my $refseq=shift @seq;
    is($seq->seq, $refseq->seq);
}
