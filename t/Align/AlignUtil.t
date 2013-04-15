# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 47 );

    use_ok( 'Bio::Align::Utilities',
        qw( aa_to_dna_aln bootstrap_replicates cat dna_to_aa_aln ) );
    use_ok('Bio::AlignIO');
    use_ok('Bio::SeqIO');
}

my $in = Bio::AlignIO->new(
    -format => 'clustalw',
    -file   => test_input_file('pep-266.aln')
);
my $pep_aln = $in->next_aln();
isa_ok( $pep_aln, 'Bio::Align::AlignI' );
$in->close();

# aa_to_dna_aln
my $seqin = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => test_input_file('cds-266.fas')
);

my %dna_seq;
while ( my $seq = $seqin->next_seq ) {
    $dna_seq{ $seq->display_id } = $seq;
}

my $dna_aln = aa_to_dna_aln( $pep_aln, \%dna_seq );

my @aa_seqs = $pep_aln->each_seq;

for my $dna_seq ( $dna_aln->each_seq ) {
    my $peptrans = $dna_seq->translate();
    my $aaseq    = shift @aa_seqs;
    is( $peptrans->seq(), $aaseq->seq() );
}

# dna_to_aa_aln
my $aa_aln = dna_to_aa_aln($dna_aln);

my @pep_seqs = $aa_aln->each_seq;

for my $dna_seq ( $dna_aln->each_seq ) {
    my $peptrans = $dna_seq->translate();
    my $aaseq    = shift @pep_seqs;
    is( $peptrans->seq, $aaseq->seq );
}

# bootstrap_replicates
my $bootstraps = bootstrap_replicates( $pep_aln, 10 );
is( scalar @$bootstraps, 10 );

# cat
my $sub_aln1 = $pep_aln->slice( 1,   100 );
my $sub_aln2 = $pep_aln->slice( 101, 200 );
my $sub_aln3 = $pep_aln->slice( 1,   200 );
my $cat_aln = cat( $sub_aln1, $sub_aln2 );
my @seq = $sub_aln3->each_seq;
for my $seq ( $cat_aln->each_seq ) {
    my $refseq = shift @seq;
    is( $seq->seq, $refseq->seq );
}
