# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 362);

    use_ok 'Bio::Seq';
    use_ok 'Bio::SeqIO';
    use_ok 'Bio::SeqFeature::Generic';
}


my ($feat, $str, $feat2, $pair, @sft); 

my $DEBUG = test_debug();

ok $feat = Bio::SeqFeature::Generic->new(
    -start => 40,
    -end => 80,
    -strand => 1,
);
is $feat->primary_tag, '';
is $feat->source_tag, '';
is $feat->display_name, '';

ok $feat = Bio::SeqFeature::Generic->new(
    -start => 40,
    -end => 80,
    -strand => 1,
    -primary => 'exon',
    -source  => 'internal',
    -display_name => 'my exon feature',
    -tag => {
        silly => 20,
        new => 1
    }
);

is $feat->start, 40, 'Start of feature location';
is $feat->end, 80, 'End of feature location';
is $feat->primary_tag, 'exon', 'Primary tag';
is $feat->source_tag, 'internal', 'Source tag';
is $feat->display_name, 'my exon feature', 'Display name';
is $feat->phase, undef, 'Undef phase by default';
is $feat->phase(1), 1, 'Phase accessor returns';
is $feat->phase, 1, 'Phase is persistent';

ok $feat->gff_string();

ok $feat2 = Bio::SeqFeature::Generic->new(
    -start => 400,
    -end => 440,
    -strand => 1,
    -primary => 'other',
    -source  => 'program_a',
        -phase => 1,
        -tag => {
            silly => 20,
            new => 1
        }
);
is $feat2->phase, 1, 'Set phase from constructor';


# Test attaching a SeqFeature::Generic to a Bio::Seq or SeqFeature::Generic
{
    # Make the parent sequence object
    my $seq = Bio::Seq->new(
        -seq        => 'aaaaggggtttt',
        -display_id => 'test',
        -alphabet   => 'dna',
    );
    
    # Make a SeqFeature
    ok my $sf1 = Bio::SeqFeature::Generic->new(
        -start  => 4,
        -end    => 9,
        -strand => 1,
    );
    
    # Add the SeqFeature to the parent
    ok $seq->add_SeqFeature($sf1);
    
    # Test that it gives the correct sequence
    is $sf1->start, 4, 'Start of first seqfeature';
    is $sf1->end, 9, 'End of first seqfeature';
    is $sf1->strand, 1, 'Strand of first seqfeature';
    ok my $sf_seq1 = $sf1->seq;
    is $sf_seq1->seq, 'aggggt', 'Sequence of first seqfeature';

    # Make a second seqfeature on the opposite strand
    ok my $sf2 = Bio::SeqFeature::Generic->new(
        -start  => 4,
        -end    => 9,
        -strand => -1,
    );
    
    # Now add the PrimarySeq to the seqfeature before adding it to the parent
    ok $sf2->attach_seq($seq->primary_seq);
    ok $seq->add_SeqFeature($sf2);
    
    # Test again that we have the correct sequence
    is $sf2->start, 4, 'Start of second seqfeature';
    is $sf2->end, 9, 'End of second seqfeature';
    is $sf2->strand, -1, 'Strand of second seqfeature';
    ok my $sf_seq2 = $sf2->seq;
    is $sf_seq2->seq, 'acccct', 'Sequence of second seqfeature';
}


# Some tests for bug #947

ok my $sfeat = Bio::SeqFeature::Generic->new(-primary => 'test');

ok $sfeat->add_sub_SeqFeature(
    Bio::SeqFeature::Generic->new(
        -start => 2,
        -end   => 4,
        -primary => 'sub1'
    ),
    'EXPAND'
);

ok $sfeat->add_sub_SeqFeature(
    Bio::SeqFeature::Generic->new(
        -start => 6,
        -end   => 8,
        -primary => 'sub2'
    ),
    'EXPAND'
);

is $sfeat->start, 2, 'sfeat start for EXPAND-ED feature (bug #947)';
is $sfeat->end, 8, 'sfeat end for EXPAND-ED feature (bug #947)';

# Some tests to see if we can set a feature to start at 0
ok $sfeat = Bio::SeqFeature::Generic->new(-start => 0, -end => 0 );

ok defined $sfeat->start;
is $sfeat->start, 0, 'Can create feature starting and ending at 0';
ok defined $sfeat->end;
is $sfeat->end, 0, 'Can create feature starting and ending at 0';


# Test for bug when Locations are not created explicitly

ok my $feat1 = Bio::SeqFeature::Generic->new(
    -start => 1,
    -end   => 15,
    -strand=> 1
);

ok $feat2 = Bio::SeqFeature::Generic->new(
    -start => 10,
    -end   => 25,
    -strand=> 1
);

ok my $overlap = $feat1->location->union($feat2->location);
is $overlap->start, 1;
is $overlap->end,   25;

ok my $intersect = $feat1->location->intersection($feat2->location);
is $intersect->start, 10;
is $intersect->end,   15;


# Now let's test spliced_seq
ok my $seqio = Bio::SeqIO->new(
    -file => test_input_file('AY095303S1.gbk'),
    -format  => 'genbank'
);
isa_ok $seqio, 'Bio::SeqIO';
ok my $geneseq = $seqio->next_seq;
isa_ok $geneseq, 'Bio::Seq';
ok my ($CDS) = grep { $_->primary_tag eq 'CDS' } $geneseq->get_SeqFeatures;
my $db;

SKIP: {
    test_skip(-tests => 5,
              -requires_modules => [qw(IO::String
                                       LWP::UserAgent
                                       HTTP::Request::Common)],
              -requires_networking => 1);
    
    use_ok 'Bio::DB::GenBank';
    
    $db = Bio::DB::GenBank->new(-verbose=> $DEBUG);
    $CDS->verbose(-1);
    my $cdsseq = $CDS->spliced_seq(-db => $db,-nosort => 1);
    
    is $cdsseq->subseq(1,76),
       'ATGCAGCCATACGCTTCCGTGAGCGGGCGATGTCTATCTAGACCAGATGCATTGCATGTGATACCGTTTGGGCGAC';
    is $cdsseq->translate->subseq(1,100),
       'MQPYASVSGRCLSRPDALHVIPFGRPLQAIAGRRFVRCFAKGGQPGDKKKLNVTDKLRLGNTPPTLDVLK'.
       'APRPTDAPSAIDDAPSTSGLGLGGGVASPR';
    # Test what happens without 
    $cdsseq = $CDS->spliced_seq(-db => $db,-nosort => 1);    
    is $cdsseq->subseq(1,76), 
       'ATGCAGCCATACGCTTCCGTGAGCGGGCGATGTCTATCTAGACCAGATGCATTGCATGTGATACCGTTTGGGCGAC';
    is $cdsseq->translate->subseq(1,100), 
       'MQPYASVSGRCLSRPDALHVIPFGRPLQAIAGRRFVRCFAKGGQPGDKKKLNVTDKLRLGNTPPTLDVLK'.
       'APRPTDAPSAIDDAPSTSGLGLGGGVASPR';
} 

ok $seqio = Bio::SeqIO->new(
    -file => test_input_file('AF032047.gbk'),
    -format  => 'genbank'
);
isa_ok $seqio, 'Bio::SeqIO';
ok $geneseq = $seqio->next_seq;
isa_ok $geneseq, 'Bio::Seq';
ok( ($CDS) = grep { $_->primary_tag eq 'CDS' } $geneseq->get_SeqFeatures );
SKIP: { 
    test_skip(-tests => 2,
              -requires_modules => [qw(IO::String
                                       LWP::UserAgent
                                       HTTP::Request::Common)],
              -requires_networking => 1);
    
    my $cdsseq = $CDS->spliced_seq( -db => $db, -nosort => 1);
    is $cdsseq->subseq(1,70),
       'ATGGCTCGCTTCGTGGTGGTAGCCCTGCTCGCGCTACTCTCTCTGTCTGGCCTGGAGGCTATCCAGCATG';
    is $cdsseq->translate->seq,
       'MARFVVVALLALLSLSGLEAIQHAPKIQVYSRHPAENGKPNFLNCYVSGFHPSDIEVDLLKNGKKIEKVE'.
       'HSDLSFSKDWSFYLLYYTEFTPNEKDEYACRVSHVTFPTPKTVKWDRTM*';
}


# Trans-spliced 

ok $seqio = Bio::SeqIO->new(
    -format => 'genbank',
    -file => test_input_file('NC_001284.gbk')
);
isa_ok $seqio, 'Bio::SeqIO';
ok my $genome = $seqio->next_seq;

for my $cds (grep { $_->primary_tag eq 'CDS' } $genome->get_SeqFeatures) {
   ok my $spliced = $cds->spliced_seq(-nosort => 1)->translate->seq;
   chop $spliced; # remove stop codon
   is $spliced, ($cds->get_tag_values('translation'))[0], 'spliced_seq translation matches expected';
}

# Spliced_seq phase 
ok my $seq = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => test_input_file('sbay_c127.fas')
)->next_seq;

ok my $sf = Bio::SeqFeature::Generic->new(
    -verbose => -1,
    -start => 263,
    -end => 721,
    -strand => 1,
    -primary => 'splicedgene'
);

ok $sf->attach_seq($seq);

my %phase_check = (
    'TTCAATGACT' => 'FNDFYSMGKS',
    'TCAATGACTT' => 'SMTSIPWVNQ',
    'GTTCAATGAC' => 'VQ*LLFHG*I',
);

for my $phase (-1..3) {
    ok my $sfseq = $sf->spliced_seq(-phase => $phase);
    ok exists $phase_check{$sfseq->subseq(1,10)};
    is $sfseq->translate->subseq(1,10), $phase_check{$sfseq->subseq(1,10)}, 'phase check';
}

# Tags
ok $sf->add_tag_value('note','n1');
ok $sf->add_tag_value('note','n2');
ok $sf->add_tag_value('comment','c1');
is_deeply [sort $sf->get_all_tags()], [sort qw(note comment)] , 'Tags found';
is_deeply [sort $sf->get_tagset_values('note')], [sort qw(n1 n2)] , 'get_tagset_values tag values found';
is_deeply [sort $sf->get_tagset_values(qw(note comment))], [sort qw(c1 n1 n2)] , 'get_tagset_values tag values for multiple tags found';
lives_ok { 
  is_deeply [sort $sf->get_tag_values('note')], [sort qw(n1 n2)] , 'get_tag_values tag values found';
} 'get_tag_values lives with tag';
lives_ok { 
  is_deeply [$sf->get_tagset_values('notag') ], [], 'get_tagset_values no tag values found';
} 'get_tagset_values lives with no tag';
throws_ok { $sf->get_tag_values('notag') } qr/tag value that does not exist/, 'get_tag_values throws with no tag';

# Circular sequence SeqFeature tests
$seq = Bio::SeqIO->new(
    -format => 'genbank',
    -file   => test_input_file('PX1CG.gb')
)->next_seq;

ok $seq->is_circular, 'Phi-X174 genome is circular';

# Retrieving the spliced sequence from any split location requires spliced_seq()

my %sf_data = (
    #       start
    'A'  => [3981, 136, 1, 1542, 'join(3981..5386,1..136)', 'ATGGTTCGTT'],
    'A*' => [4497, 136, 1, 1026, 'join(4497..5386,1..136)', 'ATGAAATCGC'],
    'B'  => [5075, 51,  1, 363,  'join(5075..5386,1..51)',  'ATGGAACAAC'],
);

ok my @split_sfs = grep {
    $_->location->isa('Bio::Location::SplitLocationI')
    } $seq->get_SeqFeatures();

is @split_sfs, 3, 'Only 3 split locations';

for my $sf (@split_sfs) {
    isa_ok $sf->location, 'Bio::Location::SplitLocationI';
    ok my ($tag) = $sf->get_tag_values('product');
    my ($start, $end, $strand, $length, $ftstring, $first_ten) = @{$sf_data{$tag}};
    
    is $sf->location->to_FTstring, $ftstring, 'Feature string';
    is $sf->spliced_seq->subseq(1,10), $first_ten, 'First ten nucleotides';
    is $sf->strand, $strand, 'Strand';
    is $sf->start, $start, 'Start';
    is $sf->end, $end, 'End';
    is $sf->length, $length, 'Expected length';
}
