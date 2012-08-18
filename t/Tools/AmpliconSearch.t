BEGIN {     
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 185);

    use_ok 'Bio::PrimarySeq';
    use_ok 'Bio::SeqFeature::Primer';
    use_ok 'Bio::Tools::AmpliconSearch';
}


my ($search, $amplicon, $seq, $seq2, $forward, $forward2, $reverse, $reverse2,
    $primer, $primer_seq, $annotated, $num_feats, $template_seq, $rna);


# Basic object

ok $search = Bio::Tools::AmpliconSearch->new(), 'Basic';
isa_ok $search, 'Bio::Tools::AmpliconSearch';


# Forward primer only

$seq = Bio::PrimarySeq->new(
   -seq => 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac',
);

$rna = Bio::PrimarySeq->new(
   -seq => 'acgAAACUUAAAGGAAUUGACGGacguacguacguGUACACACCGCCCGUacguac',
);


$forward = Bio::PrimarySeq->new(
   -seq => 'AAACTTAAAGGAATTGACGG',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -fwd_primer     => $forward,
   -attach_primers => 1,
), 'Forward primer only';
is $search->fwd_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->rev_primer, undef;
ok $template_seq = $search->template;
isa_ok $template_seq, 'Bio::Seq';
is $template_seq->seq, 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 4;
is $amplicon->end, 56;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'AAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac';
ok $primer = $amplicon->fwd_primer;
ok $primer_seq = $primer->seq;
is $primer_seq->seq, 'AAACTTAAAGGAATTGACGG';
is $primer->start, 4;
is $primer->end, 23;
is $primer->strand, 1;
is $amplicon = $search->next_amplicon, undef;


# Reverse primer only

$seq = Bio::PrimarySeq->new(
   -seq => 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac',
);
$reverse = Bio::SeqFeature::Primer->new(
   -seq => 'ACGGGCGGTGTGTAC',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -rev_primer     => $reverse,
   -attach_primers => 0,
), 'Reverse primer only';
is $search->fwd_primer, undef;
is $search->rev_primer->seq->seq, 'ACGGGCGGTGTGTAC';
is $search->template->seq, 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 1;
is $amplicon->end, 50;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGT';
is $amplicon->fwd_primer, undef;
is $amplicon->rev_primer, undef;
is $search->next_amplicon, undef;


# Forward and reverse primers, no amplicon

$seq = Bio::PrimarySeq->new(
   -seq => 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'AAACTTAAAGGAATTGACGG',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'CCCCCCCCCC',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Two primers, no match';
is $search->fwd_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->rev_primer->seq, 'CCCCCCCCCC';
is $search->template->seq, 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac';
is $search->next_amplicon, undef;


# Degenerate forward and reverse primers from file, single amplicon

ok $search = Bio::Tools::AmpliconSearch->new(
   -template    => $seq,
   -primer_file => test_input_file('forward_reverse_primers.fa'),
), 'Two degenerate primers from a file';
is $search->fwd_primer->seq, 'AAACTYAAAKGAATTGRCGG';
is $search->rev_primer->seq, 'ACGGGCGGTGTGTRC';
is $search->template->seq, 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 4;
is $amplicon->end, 50;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'AAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGT';
is $search->next_amplicon, undef;


# Same thing but with RNA template

ok $search = Bio::Tools::AmpliconSearch->new(
   -template    => $rna,
   -primer_file => test_input_file('forward_reverse_primers.fa'),
), 'RNA template';
is $search->fwd_primer->seq, 'AAACTYAAAKGAATTGRCGG';
is $search->rev_primer->seq, 'ACGGGCGGTGTGTRC';
is $search->template->seq, 'acgAAACUUAAAGGAAUUGACGGacguacguacguGUACACACCGCCCGUacguac';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 4;
is $amplicon->end, 50;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'AAACUUAAAGGAAUUGACGGacguacguacguGUACACACCGCCCGU';
is $search->next_amplicon, undef;


# Forward primer only, in sequence file

ok $search = Bio::Tools::AmpliconSearch->new(
   -primer_file => test_input_file('forward_primer.fa'),
), 'Forward primer from file';
ok $search->template($seq);
is $search->fwd_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->rev_primer, undef;
ok $template_seq = $search->template;
isa_ok $template_seq, 'Bio::Seq';
is $template_seq->seq, 'acgAAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 4;
is $amplicon->end, 56;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'AAACTTAAAGGAATTGACGGacgtacgtacgtGTACACACCGCCCGTacgtac';
is $amplicon = $search->next_amplicon, undef;


# Multiple amplicons

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtacgtacTTTTTTacgtaCCCCacgtaTTTTTTacgtacgtac',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Multiple amplicons';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 6;
is $amplicon->end, 25;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCacgtacgtacTTTTTT';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 31;
is $amplicon->end, 45;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCacgtaTTTTTT';
is $search->next_amplicon, undef;


# Amplicon on reverse strand

$seq = Bio::PrimarySeq->new(
   # Reverse-complement of previous sequence... should have same amplicons
   -seq => 'gtacgtacgtAAAAAAtacgtGGGGtacgtAAAAAAgtacgtacgtGGGGtacgt',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Amplicon on reverse strand';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 11;
is $amplicon->end, 25;
is $amplicon->strand, -1;
is $amplicon->seq->seq, 'CCCCacgtaTTTTTT';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 31;
is $amplicon->end, 50;
is $amplicon->strand, -1;
is $amplicon->seq->seq, 'CCCCacgtacgtacTTTTTT';
is $search->next_amplicon, undef;


# Overlapping amplicons (1)

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtacgtacgtacgCCCCacgtaTTTTTTacgtacgtac',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Overlapping amplicons (1)';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 25;
is $amplicon->end, 39;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCacgtaTTTTTT';
is $search->next_amplicon, undef;


# Overlapping amplicons (2)

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtacgtacgtacgTTTTTTacgtaTTTTTTacgtacgtac',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Overlapping amplicons (2)';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 6;
is $amplicon->end, 30;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCacgtacgtacgtacgTTTTTT';
is $search->next_amplicon, undef;


# Overlapping amplicons (3)

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtacgCCCCacgtacgaTTTTTTacgtaTTTTTTaaaaaaaaaa',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Overlapping amplicons (3)';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 17;
is $amplicon->end, 34;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCacgtacgaTTTTTT';
is $search->next_amplicon, undef;


# Amplicons on both strands

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtaacgtaTTTTTTacAAAAAAgtGGGGacgtaacgta',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Amplicons on both strands';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 6;
is $amplicon->end, 25;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCacgtaacgtaTTTTTT';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 28;
is $amplicon->end, 39;
is $amplicon->strand, -1;
is $amplicon->seq->seq, 'CCCCacTTTTTT';
is $search->next_amplicon, undef;


# Overlapping amplicons on both strands

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtacgtacagAAAAAATTTTTTacgtGGGGacgtacgtac',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Overlapping amplicons on both strands';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 22;
is $amplicon->end, 41;
is $amplicon->strand, -1;
is $amplicon->seq->seq, 'CCCCacgtAAAAAATTTTTT';
is $search->next_amplicon, undef;


# Annotate template

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtacgtacTTTTTTacAAAAAAgtGGGGacgtacgtac',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Annotated template';
ok $annotated = $search->annotate_template;
isa_ok $annotated, 'Bio::Seq';
$num_feats = 0;
for my $feat ( $annotated->get_SeqFeatures ) {
   $num_feats++;
   isa_ok $feat, 'Bio::SeqFeature::Amplicon';
}
is $num_feats, 2;
is $search->next_amplicon, undef;


# Update primers

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtacgtacTTTTTTacCTCTCTgtTGTGTGacgtacgtac',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CTCTCT',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'CACACA',
);
$forward2 = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse2 = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
), 'Update primers';
ok $search->fwd_primer($forward);
is $search->fwd_primer->seq, 'CTCTCT';
ok $search->rev_primer($reverse);
is $search->rev_primer->seq, 'CACACA';
ok $amplicon = $search->next_amplicon;
is $amplicon->seq->seq, 'CTCTCTgtTGTGTG';
is $search->next_amplicon, undef;
ok $search->fwd_primer($forward2);
ok $search->rev_primer($reverse2);
is $search->fwd_primer->seq, 'CCCC';
is $search->rev_primer->seq, 'AAAAAA';
ok $amplicon = $search->next_amplicon;
is $amplicon->seq->seq, 'CCCCacgtacgtacTTTTTT';
is $search->next_amplicon, undef;


# Update template

$seq = Bio::PrimarySeq->new(
   -seq => 'acgtaCCCCacgtacgtacTTTTTTa',
);
$seq2 = Bio::PrimarySeq->new(
   -seq => 'aCCCCgaTTTTTTgacgtacgtac',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Update template';
ok $search->template($seq);
is $search->template->seq, 'acgtaCCCCacgtacgtacTTTTTTa';
ok $amplicon = $search->next_amplicon;
is $amplicon->seq->seq, 'CCCCacgtacgtacTTTTTT';
is $search->next_amplicon, undef;
ok $search->template($seq2);
is $search->template->seq, 'aCCCCgaTTTTTTgacgtacgtac';
ok $amplicon = $search->next_amplicon;
is $amplicon->seq->seq, 'CCCCgaTTTTTT';
is $search->next_amplicon, undef;

