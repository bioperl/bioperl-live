# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 106 );

    use_ok('Bio::SeqIO::fastq');
    use_ok('Bio::Seq::Quality');
}

my $DEBUG = test_debug();

# round trip tests for write_fastq

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
my %example_files = ( # bug2335
    bug2335               => {
        'variant'       => 'sanger', 
        'seq'           => 'TTGGAATGTTGCAAATGGGAGGCAGTTTGAAATACTGAATAGGCCTCATC'.
                           'GAGAATGTGAAGTTTCAGTAAAGACTTGAGGAAGTTGAATGAGCTGATGA'.
                           'ATGGATATATG',
        'qual'          => '31 23 32 23 31 22 27 28 32 24 25 23 30 25 2 21 33 '.
                           '29 9 17 33 27 27 27 25 33 29 9 28 32 27 7 27 21 '.
                           '26 21 27 27 17 26 23 31 23 32 24 27 27 28 27 28 '.
                           '28 27 27 31 23 23 28 27 27 32 23 27 35 30 12 28 '.
                           '27 27 25 33 29 10 27 28 28 33 25 27 27 31 23 34 '.
                           '27 27 32 24 27 30 22 24 28 24 27 28 27 26 28 27 '.
                           '28 32 24 28 33 25 23 27 27 28 27 28 26',
        'display_id'    => 'DS6BPQV01A2G0A',
        'desc'          => undef,
        'count'         => 1
        },
    test1_sanger            => {
        'variant'       => 'sanger',
        'seq'           => 'TATTGACAATTGTAAGACCACTAAGGATTTTTGGGCGGCAGCGACTTGGA'.
                           'GCTCTTGTAAAAGCGCACTGCGTTCCTTTTCTTTATTCTTTTGATCTTGA'.
                           'GAATCTTCTAAAAATGCCGAAAAGAAATGTTGGGAAGAGAGCGTAATCAG'.
                           'TTTAGAAATGCTCTTGATGGTAGCTTTATGTTGATCCATTCTTCTGCCTC'.
                           'CTTTACGAATAAAATAGAATAAAACTCAAATGACTAATTACCTGTATTTT'.
                           'ACCTAATTTTGTGATAAAATTCAAGAAAATATGTTCGCCTTCAATAATTA'.
                           'TG',
        'qual'          => '37 37 37 37 37 37 37 37 37 37 37 40 38 40 40 37 '.
                           '37 37 39 39 40 39 39 39 39 39 37 33 33 33 33 33 '.
                           '39 39 34 29 28 28 38 39 39 39 39 39 39 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 38 '.
                           '38 29 29 29 34 38 37 37 33 33 33 33 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 37 37 37 37 37 38 38 '.
                           '38 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 '.
                           '37 37 37 37 37 37 37 37 37 34 34 34 38 38 37 37 '.
                           '37 37 37 37 37 37 37 37 40 40 40 40 38 38 38 38 '.
                           '40 40 40 38 38 38 40 40 40 40 40 40 40 40 40 40 '.
                           '38 38 38 38 38 32 25 25 25 25 30 30 31 32 32 31 '.
                           '31 31 31 31 31 31 31 31 19 19 19 19 19 22 22 31 '.
                           '31 31 31 31 31 31 31 32 32 31 32 31 31 31 31 31 '.
                           '31 25 25 25 28 28 30 30 30 30 30 31 31 32',
        'display_id'    => 'SRR005406.250',
        'desc'          => 'FB9GE3J10F6I2T length=302',
        'count'         => 250
                                },
    test2_solexa            => {
        'variant'       => 'solexa',
        'seq'           => 'GTATTATTTAATGGCATACACTCAA',
        'qual'          => '25 25 25 25 25 25 25 25 25 25 23 25 25 25 25 23 '.
                           '25 23 23 21 23 23 23 17 17',
        'display_id'    => 'SLXA-B3_649_FC8437_R1_1_1_183_714',
        'desc'          => undef,
        'count'         => 5
                                },
    test3_illumina          => {
        'variant'       => 'illumina',
        'seq'           => 'CCAAATCTTGAATTGTAGCTCCCCT',
        'qual'          => '15 19 24 15 17 24 24 24 24 24 19 24 24 21 24 24 '.
                           '20 24 24 24 24 20 18 13 19',
        'display_id'    => 'FC12044_91407_8_200_285_136',
        'desc'          => undef,
        'count'         => 25
                                },
    example                 => {  
        'variant'       => 'sanger', # TODO: guessing on the format here...
        'seq'           => 'GTTGCTTCTGGCGTGGGTGGGGGGG',
        'qual'          => '26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 '.
                           '13 22 26 18 24 18 18 18 18',
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
        'qual'          => '25 25 25 25 25 25 25 25 25 25 23 25 25 25 25 23 '.
                           '25 23 23 21 23 23 23 17 17',
        'display_id'    => 'SLXA-B3_649_FC8437_R1_1_1_183_714',
        'desc'          => undef,
        'count'         => 5
                                },
    solexa_faked            => {
        'variant'       => 'solexa',
        'seq'           => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN',
        'qual'          => '40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 '.
                           '24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 10 9 '.
                           '8 7 6 5 5 4 4 3 3 2 2 1 1',
        'display_id'    => 'slxa_0001_1_0001_01',
        'desc'          => undef,
        'count'         => 1
                                },
    tricky                  => {
        'variant'       => 'sanger', # TODO: guessing on the format here...
        'seq'           => 'TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA',
        'qual'          => '40 40 40 40 40 40 40 13 40 40 40 40 40 40 16 31 '.
                           '19 19 31 12 22 13 4 27 5 10 14 3 14 4 19 7 10 10 '.
                           '7 4',
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
        is(join(' ', @{$sample_seq->qual}),
                $example_files{$example}->{qual},
                "qual() matches $example");
    }
}


# test round-trip and conversions (single file of each type)

my @formats = qw(sanger illumina solexa);

my %conversion = (
    sanger_faked            => {
        'variant'       => 'sanger',
        'to_solexa'     => {},
        'to_illumina'   => {},
        'to_sanger'     => {},
                                },
    solexa_faked            => {
        'variant'       => 'solexa',
        'to_solexa'     => {},
        'to_illumina'   => {},
        'to_sanger'     => {},
                                },
    illumina_faked          => {
        'variant'       => 'illumina',
        'to_solexa'     => {},
        'to_illumina'   => {},
        'to_sanger'     => {},
                                },
);

#for my $example (sort keys %conversion) {
#    my $file = test_input_file('fastq', "$example.fastq");
#    my $variant = $conversion{$example}->{variant};
#    my $in = Bio::SeqIO->new(-format    => "fastq-$variant",
#                             -file      => $file,
#                             -verbose   => 2);  #strictest level
#    my $data = $in->next_dataset;
#    my $seq = Bio::Seq::Quality->new(%$data);
#    for my $format (@formats) {
#        # these are tested above already; here we are retaining the raw data,
#        # creating the Bio::Seq::Quality, writing out, then re-reading in and
#        # checking against the 
#    }
#}

# read file using format, grab first sequence via next_dataset (get raw data)
# for each 
# output to new file
# read back in using next_dataset
# check data structs using is_deeply

# test fastq errors

my %error = (
    # file name                
    error_diff_ids          => {
        variant         => 'sanger',
        exception       => qr/doesn't\smatch\sseq\sdescriptor/xms,
                                },
    error_long_qual         => {
        variant         => 'sanger',
        exception       => qr/doesn't\smatch\slength\sof\ssequence/xms,
                                },
    error_no_qual           => {
        variant         => 'sanger',
        exception       => qr/Missing\ssequence\sand\/or\squality\sdata/xms,
                                },
    error_qual_del          => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_qual_escape       => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_qual_null         => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_qual_space        => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_qual_tab          => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_qual_unit_sep     => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_qual_vtab         => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_short_qual        => {
        variant         => 'sanger',
        exception       => qr/doesn't\smatch\slength\sof\ssequence/,
                                },
    error_spaces            => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_tabs              => {
        variant         => 'sanger',
        exception       => qr/Unknown\ssymbol\swith\sASCII\svalue/xms,
                                },
    error_trunc_at_plus     => {
        variant         => 'sanger',
        exception       => qr/Missing\ssequence\sand\/or\squality\sdata/xms,
                                },
    error_trunc_at_qual     => {
        variant         => 'sanger',
        exception       => qr/Missing\ssequence\sand\/or\squality\sdata/xms,
                                },
    error_trunc_at_seq      => {
        variant         => 'sanger',
        exception       => qr/Missing\ssequence\sand\/or\squality\sdata/xms,
                                },
);

# test retrieval of raw data in a hash ref

for my $example (sort keys %error) {
    my $file = test_input_file('fastq', "$example.fastq");
    my $variant = $error{$example}->{variant};
    my $in = Bio::SeqIO->new(-format    => "fastq-$variant",
                             -file      => $file,
                             -verbose   => 2);  #strictest level
    my $ct = 0;
    throws_ok {        while (my $seq = $in->next_seq) {
            $ct++;
        } } $error{$example}->{exception}, "Exception caught for $example";
    #my $data = $in->next_dataset;
    #my $seq = Bio::Seq::Quality->new(%$data);
    #for my $format (@formats) {
    #    # these are tested above already; here we are retaining the raw data,
    #    # creating the Bio::Seq::Quality, writing out, then re-reading in and
    #    # checking against the 
    #}
}

#
#my $in_qual = Bio::SeqIO->new(
#    -file    => test_input_file( 'fastq', 'test3_illumina.fastq' ),
#    -variant => 'illumina',
#    -format  => 'fastq'
#);
#
#my $qual = $in_qual->next_dataset();
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
