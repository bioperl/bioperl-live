BEGIN {     
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 68);

    use_ok 'Bio::Tools::AmpliconSearch';
    use_ok 'Bio::PrimarySeq';
}



my ($search, $amplicon, $seq, $forward, $reverse);


# Basic object

ok $search = Bio::Tools::AmpliconSearch->new(), 'Basic';
isa_ok $search, 'Bio::Tools::AmpliconSearch';


# Forward primer only

$seq = Bio::PrimarySeq->new(
   -seq => 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'AAACTTAAAGGAATTGACGG',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -fwd_primer => $forward,
), 'Forward primer only';
is $search->fwd_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->rev_primer, undef;
is $search->template->seq, 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 4;
is $amplicon->end, 56;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa';
is $amplicon = $search->next_amplicon, undef;


# Reverse primer only

$seq = Bio::PrimarySeq->new(
   -seq => 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'ACGGGCGGTGTGTAC',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -rev_primer => $reverse,
), 'Reverse primer only';
is $search->fwd_primer, undef;
is $search->rev_primer->seq, 'ACGGGCGGTGTGTAC';
is $search->template->seq, 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 1;
is $amplicon->end, 50;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGT';
is $amplicon = $search->next_amplicon, undef;


# Forward and reverse primers, no amplicon

$reverse = Bio::PrimarySeq->new(
   -seq => 'CCCCCCCCCC',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Two primers, no match';
is $search->fwd_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->rev_primer->seq, 'CCCCCCCCCC';
is $search->template->seq, 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa';
is $amplicon = $search->next_amplicon, undef;


# Degenerate forward and reverse primers from file, single amplicon

ok $search = Bio::Tools::AmpliconSearch->new(
   -template    => $seq,
   -primer_file => test_input_file('forward_reverse_primers.fa'),
), 'Two degenerate primers from a file';
is $search->fwd_primer->seq, 'AAACTYAAAKGAATTGRCGG';
is $search->rev_primer->seq, 'ACGGGCGGTGTGTRC';
is $search->template->seq, 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 4;
is $amplicon->end, 50;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGT';
is $amplicon = $search->next_amplicon, undef;


# Multiple amplicons

$seq = Bio::PrimarySeq->new(
   -seq => 'aaaaaCCCCaaaaaaaaaaTTTTTTaaaaaCCCCaaaaaTTTTTTaaaaaaaaaa',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
);
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 6;
is $amplicon->end, 25;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCaaaaaaaaaaTTTTTT';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 31;
is $amplicon->end, 45;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCaaaaaTTTTTT';


# Amplicon on reverse strand

$seq = Bio::PrimarySeq->new(
   # Reverse-complement of previous sequence... should have same amplicons
   -seq => 'ttttttttttAAAAAAtttttGGGGtttttAAAAAAttttttttttGGGGttttt',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
);
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 31;
is $amplicon->end, 50;
is $amplicon->strand, -1;
is $amplicon->seq->seq, 'CCCCaaaaaaaaaaTTTTTT';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 11;
is $amplicon->end, 25;
is $amplicon->strand, -1;
is $amplicon->seq->seq, 'CCCCaaaaaTTTTTT';


# Multiple amplicons on muliple strands

# Overlap

