# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 10);
	
    use_ok('Bio::SeqIO');
    use_ok('Bio::Seq::PrimedSeq');
}

my ($seqio, $seq, $left, $right, $primed_seq, $left_test, $annseq, $amplicon, $returnedseq);

$seqio=Bio::SeqIO->new(-file=>test_input_file('primedseq.fa'));
$seq=$seqio->next_seq;
$left=Bio::SeqFeature::Primer->new(-seq=>'CTTTTCATTCTGACTGCAACG');
$right=Bio::SeqFeature::Primer->new(-seq=>'GGTGGTGCTAATGCGT');


ok $primed_seq = Bio::Seq::PrimedSeq->new(-seq=>$seq, -left_primer=>$left, -right_primer=>$right);
ok $left_test = $primed_seq->get_primer('left');
is $left_test,$left;
ok $annseq = $primed_seq->annotated_sequence; # should I check that this is what I think it is, or just be happy?
ok $amplicon=$primed_seq->amplicon->seq;
is uc($amplicon), uc('cttttcattctgactgcaacgGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAacgcattagcaccacc');
ok $returnedseq=$primed_seq->seq;
is $returnedseq->seq, $seq->seq;
