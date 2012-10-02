use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 37);

    use_ok 'Bio::PrimarySeq';
    use_ok 'Bio::SeqFeature::SubSeq';
}

my ($subseq, $subseq_seq, $subsubseq, $subsubseq_seq, $template);


# Basic SubSeq object

$subseq = Bio::SeqFeature::SubSeq->new();
isa_ok $subseq, 'Bio::SeqFeature::SubSeq';
isa_ok $subseq, 'Bio::SeqFeature::Generic';
is $subseq->length, undef;


# SubSeq with explicit sequence (sequence string)

ok $subseq = Bio::SeqFeature::SubSeq->new(
    -seq => 'CCCCCAAAAAGGGGGTTTTT',
);
is $subseq->length, 20;
ok $subseq_seq = $subseq->seq;
isa_ok $subseq_seq, 'Bio::PrimarySeq';
is $subseq_seq->seq, 'CCCCCAAAAAGGGGGTTTTT';


# SubSeq with explicit sequence (sequence object)

ok $subseq = Bio::SeqFeature::SubSeq->new(
    -seq => Bio::PrimarySeq->new( -seq => 'CCCCCAAAAAGGGGGTTTTT' ),
);
is $subseq->length, 20;
ok $subseq_seq = $subseq->seq;
isa_ok $subseq_seq, 'Bio::PrimarySeq';
is $subseq_seq->seq, 'CCCCCAAAAAGGGGGTTTTT';


# SubSeq with explicit sequence and coordinates

ok $subseq = Bio::SeqFeature::SubSeq->new(
    -seq    => Bio::PrimarySeq->new( -seq => 'CCCCCAAAAAGGGGGTTTTT' ),
    -start  => 11,
    -end    => 40,
    -strand => -1
);
is $subseq->length, 30;
ok $subseq_seq = $subseq->seq;
isa_ok $subseq_seq, 'Bio::PrimarySeq';
is $subseq_seq->seq, 'CCCCCAAAAAGGGGGTTTTT';


# Subseq with implicit sequence

$template = Bio::Seq->new( -seq => 'ATCGATCGATCCCCCAAAAAGGGGGTTTTTAGCTAGCTAT');

ok $subseq = Bio::SeqFeature::SubSeq->new(
    -start  => 11,
    -end    => 30,
    -strand => -1
);
is $subseq->length, 20;

ok $template->add_SeqFeature($subseq);
ok $subseq_seq = $subseq->seq;
isa_ok $subseq_seq, 'Bio::PrimarySeq';
is $subseq_seq->seq, 'AAAAACCCCCTTTTTGGGGG';


# Subseq with implicit sequence

$template = Bio::Seq->new( -seq => 'ATCGATCGATCCCCCAAAAAGGGGGTTTTTAGCTAGCTAT');

ok $subseq = Bio::SeqFeature::SubSeq->new(
    -start    => 11,
    -end      => 30,
    -strand   => -1,
    -template => $template,
);
is $subseq->length, 20;
ok $subseq_seq = $subseq->seq;
isa_ok $subseq_seq, 'Bio::PrimarySeq';
is $subseq_seq->seq, 'AAAAACCCCCTTTTTGGGGG';


# Sub SubSeq

ok $subsubseq = Bio::SeqFeature::SubSeq->new(
    -start    => 11,
    -end      => 15,
    -strand   => 1,
    -template => $subseq,
);
is $subsubseq->length, 5;
ok $subsubseq_seq = $subsubseq->seq;
isa_ok $subsubseq_seq, 'Bio::PrimarySeq';
is $subsubseq_seq->seq, 'CCCCC';


# One-liner

is( Bio::SeqFeature::SubSeq->new(-start=>11,-end=>30,-strand=>1,-template=>$template)->seq->seq, 'CCCCCAAAAAGGGGGTTTTT' );
