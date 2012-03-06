use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 21);

    use_ok 'Bio::PrimarySeq';
    use_ok 'Bio::SeqFeature::Primer';
    use_ok 'Bio::SeqFeature::Amplicon';
}

my ($amplicon, $amplicon_seq, $fwd_primer, $rev_primer, $template);


# Basic amplicon object

$amplicon = Bio::SeqFeature::Amplicon->new();
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
isa_ok $amplicon, 'Bio::SeqFeature::SubSeq';


# Amplicon with explicit sequence (sequence string)

ok $amplicon = Bio::SeqFeature::Amplicon->new(
    -seq => 'CCCCCAAAAAGGGGGTTTTT',
);
ok $amplicon_seq = $amplicon->seq;
isa_ok $amplicon_seq, 'Bio::PrimarySeq';
is $amplicon_seq->seq, 'CCCCCAAAAAGGGGGTTTTT';


# Amplicon with explicit sequence (sequence object)

$fwd_primer = Bio::SeqFeature::Primer->new(
    -start  => 1,
    -end    => 4,
    -strand => 1
);

$rev_primer = Bio::SeqFeature::Primer->new(
    -seq    => 'GATTA',
    -start  => 16,
    -end    => 20,
    -strand => -1
);

ok $amplicon = Bio::SeqFeature::Amplicon->new(
    -start      => 1,
    -end        => 20,
    -seq        => Bio::PrimarySeq->new( -seq => 'CCCCCAAAAAGGGGGTTTTT' ),
    -fwd_primer => $fwd_primer,
);

ok $amplicon->rev_primer($rev_primer);

is_deeply $amplicon->fwd_primer(), $fwd_primer;
is_deeply $amplicon->rev_primer(), $rev_primer;

ok $amplicon_seq = $amplicon->seq;
isa_ok $amplicon_seq, 'Bio::PrimarySeq';
is $amplicon_seq->seq, 'CCCCCAAAAAGGGGGTTTTT';


# Amplicon with implicit sequence

$template = Bio::Seq->new( -seq => 'ATCGATCGATCCCCCAAAAAGGGGGTTTTTAGCTAGCTAT');

ok $amplicon = Bio::SeqFeature::Amplicon->new(
    -start  => 11,
    -end    => 30,
    -strand => -1
);

ok $template->add_SeqFeature($amplicon);

ok $amplicon_seq = $amplicon->seq;
isa_ok $amplicon_seq, 'Bio::PrimarySeq';
is $amplicon_seq->seq, 'AAAAACCCCCTTTTTGGGGG';
