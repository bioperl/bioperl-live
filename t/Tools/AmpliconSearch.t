BEGIN {     
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 17);

    use_ok 'Bio::Tools::AmpliconSearch';
    use_ok 'Bio::PrimarySeq';
}





my ($search);

ok $search = Bio::Tools::AmpliconSearch->new();
isa_ok $search, 'Bio::Tools::AmpliconSearch';

my $seq = Bio::PrimarySeq->new(
   -seq => 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT'
);

my $forward = Bio::PrimarySeq->new(
   -seq => 'AAACTTAAAGGAATTGACGG'
);

my $reverse = Bio::PrimarySeq->new(
   -seq => 'GTACACACCGCCCGT'
);




ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -forward_primer => $forward,
);
is $search->get_template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
is $search->get_forward_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->get_reverse_primer, undef;




ok $search = Bio::Tools::AmpliconSearch->new(
   -template       => $seq,
   -forward_primer => $forward,
   -reverse_primer => $reverse,
);
is $search->get_template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
is $search->get_forward_primer->seq, 'AAACTTAAAGGAATTGACGG';
is $search->get_reverse_primer->seq, 'GTACACACCGCCCGT';




ok $search = Bio::Tools::AmpliconSearch->new(
   -template    => $seq,
   -primer_file => test_input_file('forward_reverse_primers.fa'),
);
is $search->get_template->seq, 'AAACTTAAAGGAATTGACGGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaGTACACACCGCCCGT';
is $search->get_forward_primer->seq, 'AAACTYAAAKGAATTGRCGG';
is $search->get_reverse_primer->seq, 'ACGGGCGGTGTGTRC';


ok $search->next_amplicon;


