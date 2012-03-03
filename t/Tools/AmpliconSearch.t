BEGIN {     
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 29);

    use_ok 'Bio::Tools::AmpliconSearch';
    use_ok 'Bio::PrimarySeq';
}



my ($search, $amplicon, $seq, $forward, $reverse);


## Basic object

#ok $search = Bio::Tools::AmpliconSearch->new(), 'Basic';
#isa_ok $search, 'Bio::Tools::AmpliconSearch';


## Forward primer only

#$seq = Bio::PrimarySeq->new(
#   -seq => 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT',
#);
#$forward = Bio::PrimarySeq->new(
#   -seq => 'AAACTTAAAGGAATTGACGG',
#);
#ok $search = Bio::Tools::AmpliconSearch->new(
#   -template       => $seq,
#   -forward_primer => $forward,
#), 'Forward primer only';
#is $search->forward_primer->seq, 'AAACTTAAAGGAATTGACGG';
#is $search->reverse_primer, undef;
#is $search->template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
#ok $amplicon = $search->next_amplicon;
#isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
#is $amplicon->start, 1;
#is $amplicon->end, 67;
#is $amplicon->strand, 1;
#is $amplicon->seq->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
#is $amplicon = $search->next_amplicon, undef;


## Forward and reverse primers, no amplicon

#$reverse = Bio::PrimarySeq->new(
#   -seq => 'GTACACACCGCCCGT',
#);
#ok $search = Bio::Tools::AmpliconSearch->new(
#   -template       => $seq,
#   -forward_primer => $forward,
#   -reverse_primer => $reverse,
#), 'Two primers, no match';
#is $search->forward_primer->seq, 'AAACTTAAAGGAATTGACGG';
#is $search->reverse_primer->seq, 'GTACACACCGCCCGT';
#is $search->template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
#is $amplicon = $search->next_amplicon, undef;


## Degenerate forward and reverse primers from file, single amplicon

#ok $search = Bio::Tools::AmpliconSearch->new(
#   -template    => $seq,
#   -primer_file => test_input_file('forward_reverse_primers.fa'),
#), 'Two degenerate primers from a file';
#is $search->forward_primer->seq, 'AAACTYAAAKGAATTGRCGG';
#is $search->reverse_primer->seq, 'ACGGGCGGTGTGTRC';
#is $search->template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
#ok $amplicon = $search->next_amplicon;
#isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
#is $amplicon->start, 1;
#is $amplicon->end, 67;
#is $amplicon->strand, 1;
#is $amplicon->seq->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
#is $amplicon = $search->next_amplicon, undef;


# Multiple amplicons

$seq = Bio::PrimarySeq->new(
   -seq => 'aaaaaCCCCCaaaaaaaaaaTTTTTaaaaaCCCCCaaaaaTTTTTaaaaaaaaaa',
);
$forward = Bio::PrimarySeq->new(
   -seq => 'CCCCC',
);
$reverse = Bio::PrimarySeq->new(
   -seq => 'AAAAA',
);
ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -forward_primer => $forward,
   -reverse_primer => $reverse,
);
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 6;
is $amplicon->end, 25;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCCaaaaaaaaaaTTTTT';
ok $amplicon = $search->next_amplicon;
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
is $amplicon->start, 31;
is $amplicon->end, 45;
is $amplicon->strand, 1;
is $amplicon->seq->seq, 'CCCCCaaaaaTTTTT';


# Amplicon on reverse strand


# Multiple amplicons on muliple strands

# Overlap

