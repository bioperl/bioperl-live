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

my $in_qual = Bio::SeqIO->new(
    -file   => test_input_file( 'fastq', 'test1_sanger.fastq' ),
    -format => 'fastq'
);
isa_ok( $in_qual, 'Bio::SeqIO' );

my $qual = $in_qual->next_seq();
isa_ok( $qual, 'Bio::Seq::Quality' );

my @quals = @{ $qual->qual() };
is( @quals, 326, 'number of qual values' );

my $qualslice = join( ',', @quals[ 25 .. 35 ] );
is( $qualslice, '30,17,17,16,16,16,16,21,18,18,20', 'qual slice' );

is( $qual->display_id, 'SRR005406.1' );
is( $qual->desc,       'FB9GE3J10GA1VT length=326' );

# Solexa, aka Illumina v1.0

# this is the test example from the MAQ script , better examples welcome!

$in_qual = Bio::SeqIO->new(
    -file   => test_input_file( 'fastq', 'test2_solexa.fastq' ),
    -format => 'fastq-solexa'
);

$qual = $in_qual->next_seq();
isa_ok( $qual, 'Bio::Seq::Quality' );

@quals = @{ $qual->qual() };
is( @quals, 25, 'number of qual values' );

$qualslice = join( ',', @quals[ 12 .. 24 ] );
is( $qualslice, '25,25,25,25,25,25,23,25,23,25,25,19,21', 'qual slice' );

is( $qual->display_id, 'SLXA-B3_649_FC8437_R1_1_1_610_79' );
is( $qual->desc,       undef );

# Illumina v1.3

$in_qual = Bio::SeqIO->new(
    -file   => test_input_file( 'fastq', 'test3_illumina.fastq' ),
    -format => 'fastq-illumina'
);

$qual = $in_qual->next_seq();
isa_ok( $qual, 'Bio::Seq::Quality' );

@quals = @{ $qual->qual() };
is( @quals, 25, 'number of qual values' );

$qualslice = join( ',', @quals[ 12 .. 22 ] );
is( $qualslice, '24,20,20,19,21,24,19,19,24,11,20', 'qual slice' );

is( $qual->display_id, 'FC12044_91407_8_200_406_24' );
is( $qual->desc,       undef );

# bug 2335

$in_qual = Bio::SeqIO->new(
    '-file'   => test_input_file( 'fastq', 'bug2335.fastq' ),
    '-format' => 'fastq-sanger'
);

$qual = $in_qual->next_seq();
isa_ok( $qual, 'Bio::Seq::Quality' );

@quals = @{ $qual->qual() };

is( @quals, 111, 'number of qual values' );

$qualslice = join( ',', @quals[ 0 .. 10 ] );
is( $qualslice, '31,23,32,23,31,22,27,28,32,24,25', 'qual slice' );

is( $qual->display_id, 'DS6BPQV01A2G0A' );
is( $qual->desc,       undef );

# raw data

$in_qual = Bio::SeqIO->new(
    -file    => test_input_file( 'fastq', 'test3_illumina.fastq' ),
    -variant => 'illumina',
    -format  => 'fastq'
);

$qual = $in_qual->next_dataset();

isa_ok( $qual, 'HASH' );
is( $qual->{-seq},         'GTTAGCTCCCACCTTAAGATGTTTA' );
is( $qual->{-raw_quality}, 'SXXTXXXXXXXXXTTSUXSSXKTMQ' );
is( $qual->{-id},          'FC12044_91407_8_200_406_24' );
is( $qual->{-desc},        '' );
is( $qual->{-descriptor},  'FC12044_91407_8_200_406_24' );
is(
    join( ',', @{ $qual->{-qual} }[ 0 .. 10 ] ),
    '19,24,24,20,24,24,24,24,24,24,24'
);

# can this be used in a constructor?

my $qualobj = Bio::Seq::Quality->new(%$qual);
is( $qualobj->seq,        'GTTAGCTCCCACCTTAAGATGTTTA' );
is( $qualobj->display_id, 'FC12044_91407_8_200_406_24' );
is( $qualobj->desc,       undef );
is(
    join( ',', @{ $qualobj->qual }[ 0 .. 10 ] ),
    '19,24,24,20,24,24,24,24,24,24,24'
);

# round trip tests for write_fastq

my %format = (
    'fastq-sanger'   => [ 'test1_sanger.fastq',   250 ],
    'fastq-solexa'   => [ 'test2_solexa.fastq',   5 ],
    'fastq-illumina' => [ 'test3_illumina.fastq', 25 ]
);

while ( my ( $variant, $data ) = each %format ) {
    my $outfile = "$variant.fastq";
    my ( $file, $total ) = @$data;
    $file = test_input_file( 'fastq', $file );
    my $in = Bio::SeqIO->new(
        -format => $variant,
        -file   => $file
    );
    my $out = Bio::SeqIO->new(
        -format => $variant,
        -file   => ">$outfile"
    );
    my ( $input_ct, $round_trip ) = ( 0, 0 );
    my $test_qual;
    while ( my $seq = $in->next_seq ) {
        $input_ct++;
        if ( $input_ct == 5 ) {
            $test_qual = $seq;
        }
        $out->write_seq($seq);
    }
    is( $input_ct, $total, $variant . " total" );
    $out->close;
    my $new_in = Bio::SeqIO->new(
        -format => $variant,
        -file   => $outfile
    );
    while ( my $seq = $new_in->next_seq ) {
        $round_trip++;
        if ( $round_trip == 5 ) {
            for my $att (qw(seq display_id desc)) {
                is( $seq->$att, $test_qual->$att, "Testing $att" );
            }
            is_deeply( $seq->qual, $test_qual->qual, "Testing qual" );
        }
    }
    is( $round_trip, $total, $variant . " total" );
}

