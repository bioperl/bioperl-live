BEGIN {     
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 109);

    use_ok 'Bio::Tools::AmpliconSearch';
    use_ok 'Bio::PrimarySeq';
    use_ok 'Bio::SeqFeature::Primer';
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
   -template   => $seq,
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
$forward = Bio::SeqFeature::Primer->new(
   -seq => 'AAACTTAAAGGAATTGACGG',
);
$reverse = Bio::SeqFeature::Primer->new(
   -seq => 'ACGGGCGGTGTGTAC',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template   => $seq,
   -rev_primer => $reverse,
), 'Reverse primer only';
is $search->fwd_primer, undef;
is $search->rev_primer->seq->seq, 'ACGGGCGGTGTGTAC';
is $search->template->seq, 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 1;
is $amplicon->end, 50;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGT';
is $search->next_amplicon, undef;


# Forward and reverse primers, no amplicon

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
is $search->template->seq, 'aaaAAACTTAAAGGAATTGACGGaaaaaaaaaaaaGTACACACCGCCCGTaaaaaa';
is $search->next_amplicon, undef;


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
is $search->next_amplicon, undef;


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
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Multiple amplicons';
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
is $search->next_amplicon, undef;


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
   -template   => $seq,
   -fwd_primer => $forward,
   -rev_primer => $reverse,
), 'Amplicon on reverse strand';
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
is $search->next_amplicon, undef;


# Overlapping amplicons (1)

$seq = Bio::PrimarySeq->new(
   -seq => 'aaaaaCCCCaaaaaaaaaaaaaaaCCCCaaaaaTTTTTTaaaaaaaaaa',
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
is $amplicon->seq->seq, 'CCCCaaaaaTTTTTT';
is $search->next_amplicon, undef;


# Overlapping amplicons (2)

$seq = Bio::PrimarySeq->new(
   -seq => 'aaaaaCCCCaaaaaaaaaaaaaaaTTTTTTaaaaaTTTTTTaaaaaaaaaa',
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
is $amplicon->seq->seq, 'CCCCaaaaaaaaaaaaaaaTTTTTT';
is $search->next_amplicon, undef;


# Overlapping amplicons (3)

$seq = Bio::PrimarySeq->new(
   -seq => 'aaaaaCCCCaaaaaaaCCCCaaaaaaaaTTTTTTaaaaaTTTTTTaaaaaaaaaa',
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
is $amplicon->seq->seq, 'CCCCaaaaaaaaTTTTTT';
is $search->next_amplicon, undef;


# Amplicons on both strands

$seq = Bio::PrimarySeq->new(
   -seq => 'aaaaaCCCCaaaaaaaaaaTTTTTTccAAAAAAttGGGGaaaaaaaaaa',
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
is $amplicon->seq->seq, 'CCCCaaaaaaaaaaTTTTTT';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 28;
is $amplicon->end, 39;
is $amplicon->strand, -1;
is $amplicon->seq->seq, 'CCCCaaTTTTTT';
is $search->next_amplicon, undef;


# Overlapping amplicons on both strands

#$seq = Bio::PrimarySeq->new(
#   -seq => 'aaaaaCCCCaaaaaaaaaaccAAAAAATTTTTTccttGGGGaaaaaaaaaa',
#);
#$forward = Bio::PrimarySeq->new(
#   -seq => 'CCCC',
#);
#$reverse = Bio::PrimarySeq->new(
#   -seq => 'AAAAAA',
#);
#ok $search = Bio::Tools::AmpliconSearch->new(
#   -template   => $seq,
#   -fwd_primer => $forward,
#   -rev_primer => $reverse,
#), 'Overlapping amplicons on both strands';
#ok $amplicon = $search->next_amplicon;
#isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
#is $amplicon->start, 6;
#is $amplicon->end, 33;
#is $amplicon->strand, 1;
#is $amplicon->seq->seq, 'CCCCaaaaaaaaaaccAAAAAATTTTTT';
#ok $amplicon = $search->next_amplicon;
#isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
#is $amplicon->start, 22;
#is $amplicon->end, 41;
#is $amplicon->strand, -1;
#is $amplicon->seq->seq, 'CCCCaaggAAAAAATTTTTT';
#is $search->next_amplicon, undef;


# attach_primers

# annotated_template

# bio::seqfeature::primer as input
