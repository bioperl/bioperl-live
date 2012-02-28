# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 140 );

    use_ok('Bio::SeqIO::fastq');
    use_ok('Bio::Seq::Quality');
}

my $DEBUG = test_debug();

# simple parsing, data conversion of fastq example files

my %example_files = (
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
    evil_wrapping           => {
        'variant'       => 'sanger', # TODO: guessing on the format here...
        'seq'           => 'AACCCGTCCCATCAAAGATTTTGGTTGGAACCCGAAAGGGTTTTGAATTC'.
                           'AAACCCCTTTCGGTTCCAACTATTCAATTGTTTAACTTTTTTTAAATTGA'.
                           'TGGTCTGTTGGACCATTTGTAATAATCCCCATCGGAATTTCTTT',
        'qual'          => '32 26 31 26 4 22 20 30 25 2 27 27 24 36 32 16 '.
                           '26 28 36 32 18 4 33 26 33 26 32 26 33 26 31 26 '.
                           '4 24 36 32 16 36 32 16 36 32 18 4 27 33 26 32 26 '.
                           '23 36 32 15 35 31 18 3 36 32 16 28 33 26 32 26 33 '.
                           '26 33 26 25 28 25 33 26 25 33 25 32 24 25 36 32 '.
                           '15 32 24 27 37 32 23 16 10 5 1 35 30 12 33 26 19 '.
                           '27 25 25 14 27 26 28 25 32 24 23 12 20 30 21 28 '.
                           '34 29 10 23 27 27 18 26 28 19 25 35 32 18 4 27 26 '.
                           '28 23 12 24 13 32 28 8 25 33 28 9',
        'display_id'    => 'SRR014849.203935',
        'desc'          => 'EIXKN4201B4HU6 length=144',
        'count'         => 3
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
        is(join(' ', map {sprintf("%.0f", $_)} @{$sample_seq->qual}),
                $example_files{$example}->{qual},
                "qual() matches $example");
        my $truncated = $sample_seq->trunc(1,10);
        is(scalar(@{$truncated->meta}), $truncated->length);
    }
}

# test round-trip and conversions (single file of each type)

my @variants = qw(sanger illumina solexa);

my %conversion = (  # check conversions, particularly solexa
    sanger_93            => {
        'variant'       => 'sanger',
        'to_solexa'     => {
          '-seq' => 'ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAN',
          '-qual' => [ (map {62} 0..31),(reverse(1..61)),1 ],
          '-raw_quality' => '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}|{zyxwvutsrqponmlkjihgfedcba`_^]\[ZYXWVUTSRQPONMLKJHGFECB@>;;',
          '-id' => 'Test',
          '-desc' => 'PHRED qualities from 93 to 0 inclusive',
          '-descriptor' => 'Test PHRED qualities from 93 to 0 inclusive'
        },
        'to_illumina'   => {
          '-seq' => 'ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAN',
          '-qual' => [ (map {62} 0..31),(reverse(0..61)) ],
          '-raw_quality' => '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}|{zyxwvutsrqponmlkjihgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@',
          '-id' => 'Test',
          '-desc' => 'PHRED qualities from 93 to 0 inclusive',
          '-descriptor' => 'Test PHRED qualities from 93 to 0 inclusive'
        },
        'to_sanger'     => {
          '-seq' => 'ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGAN',
          '-qual' => [reverse(0..93)],
          '-raw_quality' => '~}|{zyxwvutsrqponmlkjihgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@?>=<;:9876543210/.-,+*)(\'&%$#"!',
          '-id' => 'Test',
          '-desc' => 'PHRED qualities from 93 to 0 inclusive',
          '-descriptor' => 'Test PHRED qualities from 93 to 0 inclusive'
        },
    },
    solexa_faked            => {
            'variant'       => 'solexa',
        'to_solexa'     => {'-seq' => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN',
          '-qual' => [qw(40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 10 9 8 7 6 5 5 4 4 3 3 2 2 1 1)],
          '-raw_quality' => 'hgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@?>=<;',
          '-id' => 'slxa_0001_1_0001_01',
          '-desc' => '',
          '-descriptor' => 'slxa_0001_1_0001_01'
        },
        'to_illumina'   => {
          '-seq' => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN',
          '-qual' => [qw(40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 10 9 8 7 6 5 5 4 4 3 3 2 2 1 1)],
          '-namespace' => 'solexa',
          '-raw_quality' => 'hgfedcba`_^]\\[ZYXWVUTSRQPONMLKJJIHGFEEDDCCBBAA',
          '-id' => 'slxa_0001_1_0001_01',
          '-desc' => '',
          '-descriptor' => 'slxa_0001_1_0001_01'
        },
        'to_sanger'     => {
          '-seq' => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN',
          '-qual' => [qw(40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 10 9 8 7 6 5 5 4 4 3 3 2 2 1 1)],
          '-namespace' => 'solexa',
          '-raw_quality' => 'IHGFEDCBA@?>=<;:9876543210/.-,++*)(\'&&%%$$##""',
          '-id' => 'slxa_0001_1_0001_01',
          '-desc' => '',
          '-descriptor' => 'slxa_0001_1_0001_01'
        },
    },
    illumina_faked          => {
        'variant'       => 'illumina',
        'to_solexa'     => {
          '-seq' => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN',
          '-qual' => [reverse(1..40), 1],  # round trip from solexa is lossy
          '-namespace' => 'illumina',
          '-raw_quality' => 'hgfedcba`_^]\[ZYXWVUTSRQPONMLKJHGFECB@>;;',
          '-id' => 'Test',
          '-desc' => 'PHRED qualities from 40 to 0 inclusive',
          '-descriptor' => 'Test PHRED qualities from 40 to 0 inclusive'
        },
        'to_illumina'   => {
          '-seq' => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN',
          '-qual' => [reverse(0..40)],
          '-raw_quality' => 'hgfedcba`_^]\\[ZYXWVUTSRQPONMLKJIHGFEDCBA@',
          '-id' => 'Test',
          '-desc' => 'PHRED qualities from 40 to 0 inclusive',
          '-descriptor' => 'Test PHRED qualities from 40 to 0 inclusive'
        },
        'to_sanger'     => {
          '-seq' => 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN',
          '-qual' => [reverse(0..40)],
          '-raw_quality' => 'IHGFEDCBA@?>=<;:9876543210/.-,+*)(\'&%$#"!',
          '-id' => 'Test',
          '-desc' => 'PHRED qualities from 40 to 0 inclusive',
          '-descriptor' => 'Test PHRED qualities from 40 to 0 inclusive'
        }
    },
);

for my $example (sort keys %conversion) {
    my $file = test_input_file('fastq', "$example.fastq");
    my $variant = $conversion{$example}->{variant};
    my $in = Bio::SeqIO->new(-format    => "fastq-$variant",
                             -file      => $file,
                             -verbose   => 2);  #strictest level
    # this both tests the next_dataset method and helps check roundtripping
    my $seq = $in->next_seq;
    for my $newvar (@variants) {
        next unless exists $conversion{$example}->{"to_$newvar"};
        my $outfile = test_output_file();
        Bio::SeqIO->new(-format   => "fastq-$newvar",
                        -file     => ">$outfile",
                        -verbose  => -1)->write_seq($seq);
        my $newdata = Bio::SeqIO->new(-format => "fastq-$newvar",
                                    -file     => $outfile)->next_dataset;
        # round for simple comparison, get around floating pt comparison probs

        if ($newvar eq 'solexa') {
            $newdata->{-qual} = [map {sprintf("%.0f",$_)} @{$newdata->{-qual}}];
        }

        #print Dumper($newdata) if $variant eq 'sanger' && $newvar eq 'illumina';

        $conversion{$example}->{"to_$newvar"}->{'-namespace'} = $newvar;
        is_deeply($newdata, $conversion{$example}->{"to_$newvar"}, "Conversion from $variant to $newvar");
    }
}

# test fastq exception handling

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
    error_trunc_in_title    => {
        variant         => 'sanger',
        exception       => qr/Missing\ssequence\sand\/or\squality\sdata/xms,
                                },
    error_trunc_in_seq      => {
        variant         => 'sanger',
        exception       => qr/Missing\ssequence\sand\/or\squality\sdata/xms,
                                },
    error_trunc_in_plus     => {
        variant         => 'sanger',
        exception       => qr/doesn't\smatch\sseq\s descriptor/xms,
                                },
    error_trunc_in_qual     => {
        variant         => 'sanger',
        exception       => qr/doesn't\smatch\slength\sof\ssequence/xms,
                                },
);

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
}

# fastq

my $in = Bio::SeqIO->new(-format    => 'fastq',
                         -file      => test_input_file('fastq', 'zero_qual.fastq'),
                         -verbose   => 2);  # strictest level

lives_and {my $seq = $in->next_seq;
           is($seq->seq, 'G');} 'edge case; single 0 in quality fails';


