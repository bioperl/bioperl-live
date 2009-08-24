# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 52 );

    use_ok('Bio::SeqIO::fastq');
    use_ok('Bio::Seq::Quality');
}

my $DEBUG = test_debug();

# original FASTQ (Sanger); data is from NCBI SRA database, which has
# all data converted over to Sanger version of FASTQ

#my $in_qual = Bio::SeqIO->new(
#    -file   => test_input_file( 'fastq', 'test1_sanger.fastq' ),
#    -format => 'fastq'
#);
#isa_ok( $in_qual, 'Bio::SeqIO' );
#
#my $qual = $in_qual->next_seq();
#isa_ok( $qual, 'Bio::Seq::Quality' );
#
#my @quals = @{ $qual->qual() };
#is( @quals, 326, 'number of qual values' );
#
#my $qualslice = join( ',', @quals[ 25 .. 35 ] );
#is( $qualslice, '30,17,17,16,16,16,16,21,18,18,20', 'qual slice' );
#
#is( $qual->display_id, 'SRR005406.1' );
#is( $qual->desc,       'FB9GE3J10GA1VT length=326' );
#
## Solexa, aka Illumina v1.0
#
## this is the test example from the MAQ script , better examples welcome!
#
#$in_qual = Bio::SeqIO->new(
#    -file   => test_input_file( 'fastq', 'test2_solexa.fastq' ),
#    -format => 'fastq-solexa'
#);
#
#$qual = $in_qual->next_seq();
#isa_ok( $qual, 'Bio::Seq::Quality' );
#
#@quals = @{ $qual->qual() };
#is( @quals, 25, 'number of qual values' );
#
#$qualslice = join( ',', @quals[ 12 .. 24 ] );
#is( $qualslice, '25,25,25,25,25,25,23,25,23,25,25,19,21', 'qual slice' );
#
#is( $qual->display_id, 'SLXA-B3_649_FC8437_R1_1_1_610_79' );
#is( $qual->desc,       undef );
#
## Illumina v1.3
#
#$in_qual = Bio::SeqIO->new(
#    -file   => test_input_file( 'fastq', 'test3_illumina.fastq' ),
#    -format => 'fastq-illumina'
#);
#
#$qual = $in_qual->next_seq();
#isa_ok( $qual, 'Bio::Seq::Quality' );
#
#@quals = @{ $qual->qual() };
#is( @quals, 25, 'number of qual values' );
#
#$qualslice = join( ',', @quals[ 12 .. 22 ] );
#is( $qualslice, '24,20,20,19,21,24,19,19,24,11,20', 'qual slice' );
#
#is( $qual->display_id, 'FC12044_91407_8_200_406_24' );
#is( $qual->desc,       undef );
#
## bug 2335
#
#$in_qual = Bio::SeqIO->new(
#    '-file'   => test_input_file( 'fastq', 'bug2335.fastq' ),
#    '-format' => 'fastq-sanger'
#);
#
#$qual = $in_qual->next_seq();
#isa_ok( $qual, 'Bio::Seq::Quality' );
#
#@quals = @{ $qual->qual() };
#
#is( @quals, 111, 'number of qual values' );
#
#$qualslice = join( ',', @quals[ 0 .. 10 ] );
#is( $qualslice, '31,23,32,23,31,22,27,28,32,24,25', 'qual slice' );
#
#is( $qual->display_id, 'DS6BPQV01A2G0A' );
#is( $qual->desc,       undef );
#
## raw data
#
#$in_qual = Bio::SeqIO->new(
#    -file    => test_input_file( 'fastq', 'test3_illumina.fastq' ),
#    -variant => 'illumina',
#    -format  => 'fastq'
#);
#
#$qual = $in_qual->next_dataset();
#
#isa_ok( $qual, 'HASH' );
#is( $qual->{-seq},         'GTTAGCTCCCACCTTAAGATGTTTA' );
#is( $qual->{-raw_quality}, 'SXXTXXXXXXXXXTTSUXSSXKTMQ' );
#is( $qual->{-id},          'FC12044_91407_8_200_406_24' );
#is( $qual->{-desc},        '' );
#is( $qual->{-descriptor},  'FC12044_91407_8_200_406_24' );
#is(
#    join( ',', @{ $qual->{-qual} }[ 0 .. 10 ] ),
#    '19,24,24,20,24,24,24,24,24,24,24'
#);
#
## can this be used in a constructor?
#
#my $qualobj = Bio::Seq::Quality->new(%$qual);
#is( $qualobj->seq,        'GTTAGCTCCCACCTTAAGATGTTTA' );
#is( $qualobj->display_id, 'FC12044_91407_8_200_406_24' );
#is( $qualobj->desc,       undef );
#is(
#    join( ',', @{ $qualobj->qual }[ 0 .. 10 ] ),
#    '19,24,24,20,24,24,24,24,24,24,24'
#);
#
## round trip tests for write_fastq
#
#my %format = (
#    'fastq-sanger'   => [ 'test1_sanger.fastq',   250 ],
#    'fastq-solexa'   => [ 'test2_solexa.fastq',   5 ],
#    'fastq-illumina' => [ 'test3_illumina.fastq', 25 ]
#);
#
#while ( my ( $variant, $data ) = each %format ) {
#    my $outfile = test_output_file();
#    my ( $file, $total ) = @$data;
#    $file = test_input_file( 'fastq', $file );
#    my $in = Bio::SeqIO->new(
#        -format => $variant,
#        -file   => $file
#    );
#    my $out = Bio::SeqIO->new(
#        -format => $variant,
#        -file   => ">$outfile"
#    );
#    my ( $input_ct, $round_trip ) = ( 0, 0 );
#    my $test_qual;
#    while ( my $seq = $in->next_seq ) {
#        $input_ct++;
#        if ( $input_ct == 5 ) {
#            $test_qual = $seq;
#        }
#        $out->write_seq($seq);
#    }
#    is( $input_ct, $total, $variant . " total" );
#    $out->close;
#    my $new_in = Bio::SeqIO->new(
#        -format => $variant,
#        -file   => $outfile
#    );
#    while ( my $seq = $new_in->next_seq ) {
#        $round_trip++;
#        if ( $round_trip == 5 ) {
#            for my $att (qw(seq display_id desc)) {
#                is( $seq->$att, $test_qual->$att, "Testing $att" );
#            }
#            is_deeply( $seq->qual, $test_qual->qual, "Testing qual" );
#        }
#    }
#    is( $round_trip, $total, $variant . " total" );
#}

# test simple parsing of fastq example files
my %example_files = (
    example                 => {
        'variant'       => 'sanger',  # not sure about this one
        'seq'           => 'GTTGCTTCTGGCGTGGGTGGGGGGG',
        'qual'          => '',
        'display_id'    => 'EAS54_6_R1_2_1_443_348',
        'desc'          => undef,
        'count'         => 3
                                },
    illumina_faked          => {
        'variant'       => 'illumina',
        'seq'           => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN',
        'qual'          => '40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 '.
                           '21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0',
        'display_id'    => 'Test',
        'desc'          => 'PHRED qualities from 40 to 0 inclusive',
        'count'         => 1
                                },
    sanger_93               => {
        'variant'       => 'sanger',
        'seq'           => 'ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAC'.
                           'TGACTGACTGACTGACTGACTGACTGACTGACTGAN',
        'qual'          => '93 92 91 90 89 88 87 86 85 84 83 82 81 80 79 78 77 76 75 '.
                           '74 73 72 71 70 69 68 67 66 65 64 63 62 61 60 59 58 57 56 '.
                           '55 54 53 52 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 '.
                           '36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 '.
                           '17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0',
        'display_id'    => 'Test',
        'desc'          => 'PHRED qualities from 93 to 0 inclusive',
        'count'         => 1
                                },
    sanger_faked            => {
        'variant'       => 'sanger',
        'seq'           => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN',
        'qual'          => '40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 '.
                            '21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0',
        'display_id'    => 'Test',
        'desc'          => 'PHRED qualities from 40 to 0 inclusive',
        'count'         => 1
                                },
    solexa_example          => {
        'variant'       => 'solexa',
        'seq'           => 'GTATTATTTAATGGCATACACTCAA',
        'qual'          => '',
        'display_id'    => 'SLXA-B3_649_FC8437_R1_1_1_183_714',
        'desc'          => undef,
        'count'         => 5
                                },
    solexa_faked            => {
        'variant'       => 'solexa',
        'seq'           => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN',
        'qual'          => '',
        'display_id'    => 'slxa_0001_1_0001_01',
        'desc'          => undef,
        'count'         => 1
                                },
    tricky                  => {
        'variant'       => 'solexa', # not sure about this one
        'seq'           => 'TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA',
        'qual'          => '',
        'display_id'    => '071113_EAS56_0053:1:3:990:501',
        'desc'          => undef,
        'count'         => 4
                                },
);

for my $example (sort keys %example_files) {
    my $file = test_input_file('fastq', "$example.fastq");
    my $variant = $example_files{$example}->{variant};
    my $in = Bio::SeqIO->new(-format    => "fastq-$variant",
                             -file      => $file,
                             -verbose   => 2);  #strictest level
    my $ct = 0;
    my $sample_seq;
    eval {
        while (my $seq = $in->next_seq) {
            $ct++;
            $sample_seq = $seq; # always grab the last seq
        }
    };
    ok(!$@, "$example parses");
    is($ct, $example_files{$example}->{count}, "correct num. seqs in $example");
    ok(defined($sample_seq), 'sample sequence obtained');
    if ($sample_seq) {
        isa_ok($sample_seq, 'Bio::Seq::Quality');
        for my $method (qw(seq desc display_id)) {
            is($sample_seq->$method,
               $example_files{$example}->{$method},
               "$method() matches $example");
        }
        is(join(' ',@{$sample_seq->qual}),
                $example_files{$example}->{qual},
                "qual() matches $file");
    }
}


# test conversions (single files of each type)

my %conversions = (
    illumina_faked          => {'variant'       => 'solexa',
                                'seq'           => '',
                                'qual'          => '',
                                'display_id'    => '',
                                'desc'          => ''},
    sanger_faked            => {'variant'       => 'solexa',
                                'seq'           => '',
                                'qual'          => '',
                                'display_id'    => '',
                                'desc'          => ''},
    solexa_faked            => {'variant'       => 'solexa',
                                'seq'           => '',
                                'qual'          => '',
                                'display_id'    => '',
                                'desc'          => ''},
);

# test fastq errors/warnings

my %error = (
    # file name                exception
    error_diff_ids          => qr//,
    error_long_qual         => qr//,
    error_no_qual           => qr//,
    error_qual_del          => qr//,
    error_qual_escape       => qr//,
    error_qual_null         => qr//,
    error_qual_space        => qr//,
    error_qual_tab          => qr//,
    error_qual_unit_sep     => qr//,
    error_qual_vtab         => qr//,
    error_short_qual        => qr//,
    error_spaces            => qr//,
    error_tabs              => qr//,
    error_trunc_at_plus     => qr//,
    error_trunc_at_qual     => qr//,
    error_trunc_at_seq      => qr//
);

