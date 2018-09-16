# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 65);

    use_ok('Bio::SeqIO');
    use_ok('Bio::Seq::PrimedSeq');
}

my ($seqio, $seq, $left, $right, $primed_seq, $left_test, $right_test, $annseq, $amplicon, $returnedseq);

$seqio = Bio::SeqIO->new(-file => test_input_file('primedseq.fa'));
$seq   = $seqio->next_seq;

my $expected_amplicon_seq = 'cttttcattctgactgcaacgGGCAATATGTCTCTGTGTGGATTAAAAA'.
'AAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAA'.
'TACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAacgcattagca'.
'ccacc';

my $expected_amplicon_seq2 = 'cttttcattctgactgcaacgTGTCTCTGTGTGGATTAAAAA'.
'AAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAA'.
'TACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAacgcattagca'.
'ccacc';

# Prime with Bio::PrimarySeqI objects and have the primer positions calculated
$left  = Bio::PrimarySeq->new(-id => 123, -seq => 'CTTTTCATTCTGACTGCAACG');
$right = Bio::Seq->new(-seq => 'GGTGGTGCTAATGCGT');
ok $primed_seq = Bio::Seq::PrimedSeq->new(
    -seq          => $seq,
    -left_primer  => $left,
    -right_primer => $right,
), 'Priming the target with sequence objects';

is $primed_seq->isa('Bio::SeqFeature::Generic'), 1;

ok $annseq = $primed_seq->annotated_sequence; # should I check that this is what I think it is, or just be happy?

ok $amplicon = $primed_seq->amplicon;
is $amplicon->seq, $expected_amplicon_seq;
is $amplicon->id, 'Amplicon_from_Test1';
ok $returnedseq = $primed_seq->seq;
is $returnedseq->seq, $seq->seq;
ok $left_test = $primed_seq->get_primer('-left_primer');
isa_ok $left_test, 'Bio::SeqFeature::Primer';
is $left_test->seq->seq, 'CTTTTCATTCTGACTGCAACG';
ok $right_test = $primed_seq->get_primer('-r');
isa_ok $right_test, 'Bio::SeqFeature::Primer';
is $right_test->seq->seq, 'GGTGGTGCTAATGCGT';
ok( ($left_test, $right_test) = $primed_seq->get_primer() );
is $left_test->seq->seq, 'CTTTTCATTCTGACTGCAACG';
is $right_test->seq->seq, 'GGTGGTGCTAATGCGT';
is $left_test->strand, 1;
is $left_test->start, 3;
is $left_test->end, 23;
is $right_test->strand, -1;
is $right_test->start, 195;
is $right_test->end, 210;


# Prime the sequence with Bio::SeqFeature::Primer objects
$left  = Bio::SeqFeature::Primer->new(-id => 123, -seq => 'CTTTTCATTCTGACTGCAACG');
$right = Bio::SeqFeature::Primer->new(-seq => 'GGTGGTGCTAATGCGT');
ok $primed_seq = Bio::Seq::PrimedSeq->new(
    -seq          => $seq,
    -left_primer  => $left,
    -right_primer => $right,
), 'Priming the target with primer objects';
ok $annseq = $primed_seq->annotated_sequence;
ok $amplicon = $primed_seq->amplicon;
is $amplicon->seq, $expected_amplicon_seq;
is $amplicon->id, 'Amplicon_from_Test1';
ok $returnedseq = $primed_seq->seq;
is $returnedseq->seq, $seq->seq;
ok $left_test = $primed_seq->get_primer('left');
is_deeply $left_test, $left;
ok $right_test = $primed_seq->get_primer('r');
is_deeply $right_test, $right;
ok( ($left_test, $right_test) = $primed_seq->get_primer('-both') );
is_deeply $left_test, $left;
is_deeply $right_test, $right;
is $left_test->strand, 1;
is $left_test->start, 3;
is $left_test->end, 23;
is $right_test->strand, -1;
is $right_test->start, 195;
is $right_test->end, 210;


# Prime the sequence with Bio::SeqFeature::Primer objects
$left  = Bio::SeqFeature::Primer->new(
    -id => 123,
    -seq => 'CTTTTCATTCTGACTGCAACG',
    -start => 10,
    -end => 30,
    -strand => 1,
);
$right = Bio::SeqFeature::Primer->new(
    -seq => 'GGTGGTGCTAATGCGT',
    -start => 195,
    -end => 210,
    -strand => -1,
);
ok $primed_seq = Bio::Seq::PrimedSeq->new(
    -seq          => $seq,
    -left_primer  => $left,
    -right_primer => $right,
), 'Priming the target with located primer objects';
ok $annseq = $primed_seq->annotated_sequence;
ok $amplicon = $primed_seq->amplicon;
is $amplicon->seq, $expected_amplicon_seq2;
is $amplicon->id, 'Amplicon_from_Test1';
ok $returnedseq = $primed_seq->seq;
is $returnedseq->seq, $seq->seq;
ok $left_test = $primed_seq->get_primer('left');
is_deeply $left_test, $left;
ok $right_test = $primed_seq->get_primer('r');
is_deeply $right_test, $right;
ok( ($left_test, $right_test) = $primed_seq->get_primer('-both') );
is_deeply $left_test, $left;
is_deeply $right_test, $right;
is $left_test->strand, 1;
is $left_test->start, 10;
is $left_test->end, 30;
is $right_test->strand, -1;
is $right_test->start, 195;
is $right_test->end, 210;
