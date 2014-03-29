# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_blast.t 14995 2008-11-16 06:20:00Z cjfields $

use strict;
use warnings;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 1389);

    use_ok('Bio::SearchIO');
}

SKIP: {
    test_skip(-tests => 4, -requires_module => 'Path::Class');


    my $file = Path::Class::file(test_input_file('ecolitst.bls'));
    my $f    = sub { my ($file) = @_; Bio::SearchIO->new( -file   => $file, -format => 'blast') };

    lives_ok(sub { $f->($file) } , 'Bio::SearchIO->new can handle a Path::Class object');
    isa_ok($f->($file), 'Bio::Root::IO');

    $file = Path::Class::dir(File::Spec->catfile(qw/t data/))->file('ecolitst.bls');

    lives_ok(sub { $f->($file) } , 'Bio::SearchIO->new can handle a Path::Class object');
    isa_ok($f->($file), 'Bio::Root::IO');
}

my ( $searchio, $result, $iter, $hit, $hsp );

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('ecolitst.bls')
);

$result = $searchio->next_result;

like($result->algorithm_reference,
    qr/Gapped BLAST and PSI-BLAST: a new generation of protein database search/
);

is( $result->database_name, 'ecoli.aa', 'database_name()' );
is( $result->database_entries, 4289 );
is( $result->database_letters, 1358990 );

is( $result->algorithm, 'BLASTP' );
like( $result->algorithm_version, qr/^2\.1\.3/ );
like( $result->query_name,
qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/
);
is( $result->query_accession,                 'AAC73113.1' );
is( $result->query_gi,                        1786183 );
is( $result->query_length,                    820 );
is( $result->get_statistic('kappa'),          '0.135' );
is( $result->get_statistic('kappa_gapped'),   '0.0410' );
is( $result->get_statistic('lambda'),         '0.319' );
is( $result->get_statistic('lambda_gapped'),  '0.267' );
is( $result->get_statistic('entropy'),        '0.383' );
is( $result->get_statistic('entropy_gapped'), '0.140' );

is( $result->get_statistic('dbletters'),           1358990 );
is( $result->get_statistic('dbentries'),           4289 );
is( $result->get_statistic('effective_hsplength'), 47 );
is( $result->get_statistic('effectivespace'),      894675611 );
is( $result->get_parameter('matrix'),              'BLOSUM62' );
is( $result->get_parameter('gapopen'),             11 );
is( $result->get_parameter('gapext'),              1 );
is( $result->get_statistic('S2'),                  '92' );
is( $result->get_statistic('S2_bits'),             '40.0' );
float_is( $result->get_parameter('expect'), '1.0e-03' );
is( $result->get_statistic('num_extensions'),     '82424' );
is( $result->get_statistic('querylength'),        773 );
is( $result->get_statistic('effectivedblength'),  1157407 );
is( $result->get_statistic('effectivespaceused'), 894675611 );

my @valid = (
    [ 'gb|AAC73113.1|', 820, 'AAC73113', '0',     1567, 4058 ],
    [ 'gb|AAC76922.1|', 810, 'AAC76922', '1e-91', 332,  850 ],
    [ 'gb|AAC76994.1|', 449, 'AAC76994', '3e-47', 184,  467 ]
);
my $count = 0;
while ( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->bits,      shift @$d );
    is( $hit->raw_score, shift @$d );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    1 );
            is( $hsp->query->end,      820 );
            is( $hsp->hit->start,      1 );
            is( $hsp->hit->end,        820 );
            is( $hsp->length('total'), 820 );
            is( $hsp->start('hit'),    $hsp->hit->start );
            is( $hsp->end('query'),    $hsp->query->end );
            is( $hsp->strand('sbjct'), $hsp->subject->strand );  # alias for hit
            float_is( $hsp->evalue, 0.0 );
            is( $hsp->score, 4058 );
            is( $hsp->bits,  1567 );
            is( sprintf( "%.2f", $hsp->percent_identity ),        98.29 );
            is( sprintf( "%.4f", $hsp->frac_identical('query') ), 0.9829 );
            is( sprintf( "%.4f", $hsp->frac_identical('hit') ),   0.9829 );
            is( $hsp->gaps, 0 );
            is( $hsp->n,    1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( @valid, 0 );

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('ecolitst.wublastp')
);

$result = $searchio->next_result;

like($result->algorithm_reference,
     qr/Gish, W. \(1996-2000\)/);

is( $result->database_name,    'ecoli.aa' );
is( $result->database_letters, 1358990 );
is( $result->database_entries, 4289 );
is( $result->algorithm,        'BLASTP' );
like( $result->algorithm_version, qr/^2\.0MP\-WashU/ );
like( $result->query_name,
qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/
);
is( $result->query_accession, 'AAC73113.1' );

is( $result->query_length,                          820 );
is( $result->query_gi,                              1786183 );
is( $result->get_statistic('kappa'),                0.136 );
is( $result->get_statistic('lambda'),               0.319 );
is( $result->get_statistic('entropy'),              0.384 );
is( $result->get_statistic('dbletters'),            1358990 );
is( $result->get_statistic('dbentries'),            4289 );
is( $result->get_parameter('matrix'),               'BLOSUM62' );
is( $result->get_statistic('Frame+0_lambda_used'),  '0.319' );
is( $result->get_statistic('Frame+0_kappa_used'),   '0.136' );
is( $result->get_statistic('Frame+0_entropy_used'), '0.384' );

is( $result->get_statistic('Frame+0_lambda_computed'),  '0.319' );
is( $result->get_statistic('Frame+0_kappa_computed'),   '0.136' );
is( $result->get_statistic('Frame+0_entropy_computed'), '0.384' );

is( $result->get_statistic('Frame+0_lambda_gapped'),  '0.244' );
is( $result->get_statistic('Frame+0_kappa_gapped'),   '0.0300' );
is( $result->get_statistic('Frame+0_entropy_gapped'), '0.180' );

@valid = (
    [ 'gb|AAC73113.1|', 820, 'AAC73113', '0',       4141 ],
    [ 'gb|AAC76922.1|', 810, 'AAC76922', '3.1e-86', 844 ],
    [ 'gb|AAC76994.1|', 449, 'AAC76994', '2.8e-47', 483 ]
);
$count = 0;
while ( $hit = $result->next_hit ) {
    my $d = shift @valid;

    if ( $count == 1 ) {

        # Test HSP contig data returned by SearchUtils::tile_hsps()
        # Second hit has two hsps that overlap.

        # compare with the contig made by hand for these two contigs
        # in t/data/contig-by-hand.wublastp
        # (in this made-up file, the hsps from ecolitst.wublastp
        #  were aligned and contiged, and Length, Identities, Positives
        #  were counted, by a human (maj) )

        my $hand_hit = Bio::SearchIO->new(
            -format => 'blast',
            -file   => test_input_file('contig-by-hand.wublastp')
        )->next_result->next_hit;
        my $hand_hsp     = $hand_hit->next_hsp;
        my @hand_qrng    = $hand_hsp->range('query');
        my @hand_srng    = $hand_hsp->range('hit');
        my @hand_matches = $hand_hit->matches;

        my ( $qcontigs, $scontigs ) = Bio::Search::SearchUtils::tile_hsps($hit);

        # Query contigs
        is( $qcontigs->[0]->{'start'}, $hand_qrng[0] );
        is( $qcontigs->[0]->{'stop'},  $hand_qrng[1] );
        is( $qcontigs->[0]->{'iden'},  $hand_matches[0] );
        is( $qcontigs->[0]->{'cons'},  $hand_matches[1] );

        # Subject contigs
        is( $scontigs->[0]->{'start'}, $hand_srng[0] );
        is( $scontigs->[0]->{'stop'},  $hand_srng[1] );
        is( $scontigs->[0]->{'iden'},  $hand_matches[0] );
        is( $scontigs->[0]->{'cons'},  $hand_matches[1] );
    }

    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->raw_score, shift @$d );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    1 );
            is( $hsp->query->end,      820 );
            is( $hsp->hit->start,      1 );
            is( $hsp->hit->end,        820 );
            is( $hsp->length('total'), 820 );

            float_is( $hsp->evalue, 0.0 );
            float_is( $hsp->pvalue, '0.0' );
            is( $hsp->score,                   4141 );
            is( $hsp->bits,                    1462.8 );
            is( $hsp->percent_identity,        100 );
            is( $hsp->frac_identical('query'), 1.00 );
            is( $hsp->frac_identical('hit'),   1.00 );
            is( $hsp->gaps,                    0 );
            is( $hsp->n,                       1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( @valid, 0 );

# test that add hit really works properly for BLAST objects
# bug 1611
my @hits = $result->hits;
$result->add_hit( $hits[0] );
is( $result->num_hits, @hits + 1 );

# test WU-BLAST -noseqs option
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('ecolitst.noseqs.wublastp')
);

$result = $searchio->next_result;
is(
    $result->algorithm_reference, 'Gish, W. (1996-2004) http://blast.wustl.edu
'
);
is( $result->database_name,    'ecoli.aa' );
is( $result->database_letters, 1358990 );
is( $result->database_entries, 4289 );
is( $result->algorithm,        'BLASTP' );
like( $result->algorithm_version, qr/^2\.0MP\-WashU/ );
like( $result->query_name,
qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/
);
is( $result->query_accession, 'AAC73113.1' );
is( $result->query_gi,        1786183 );

is( $result->query_length,                          820 );
is( $result->get_statistic('kappa'),                0.135 );
is( $result->get_statistic('lambda'),               0.319 );
is( $result->get_statistic('entropy'),              0.384 );
is( $result->get_statistic('dbletters'),            1358990 );
is( $result->get_statistic('dbentries'),            4289 );
is( $result->get_parameter('matrix'),               'BLOSUM62' );
is( $result->get_statistic('Frame+0_lambda_used'),  '0.319' );
is( $result->get_statistic('Frame+0_kappa_used'),   '0.135' );
is( $result->get_statistic('Frame+0_entropy_used'), '0.384' );

is( $result->get_statistic('Frame+0_lambda_computed'),  '0.319' );
is( $result->get_statistic('Frame+0_kappa_computed'),   '0.135' );
is( $result->get_statistic('Frame+0_entropy_computed'), '0.384' );

is( $result->get_statistic('Frame+0_lambda_gapped'),  '0.244' );
is( $result->get_statistic('Frame+0_kappa_gapped'),   '0.0300' );
is( $result->get_statistic('Frame+0_entropy_gapped'), '0.180' );

@valid = (
    [ 'gb|AAC73113.1|', 820, 'AAC73113', '0',       4141 ],
    [ 'gb|AAC76922.1|', 810, 'AAC76922', '6.6e-93', 907 ],
    [ 'gb|AAC76994.1|', 449, 'AAC76994', '2.8e-47', 483 ]
);
$count = 0;
while ( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->raw_score, shift @$d );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    1 );
            is( $hsp->query->end,      820 );
            is( $hsp->hit->start,      1 );
            is( $hsp->hit->end,        820 );
            is( $hsp->length('total'), 820 );

            float_is( $hsp->evalue, 0. );
            float_is( $hsp->pvalue, '0.' );
            is( $hsp->score,                   4141 );
            is( $hsp->bits,                    1462.8 );
            is( $hsp->percent_identity,        100 );
            is( $hsp->frac_identical('query'), 1.00 );
            is( $hsp->frac_identical('hit'),   1.00 );
            is( $hsp->gaps,                    0 );
            is( $hsp->n,                       1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( @valid, 0 );

# test tblastx
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('HUMBETGLOA.tblastx')
);

$result = $searchio->next_result;
like($result->algorithm_reference,qr/Gapped BLAST and PSI-BLAST/);
is( $result->database_name,    'ecoli.nt' );
is( $result->database_letters, 4662239 );
is( $result->database_entries, 400 );
is( $result->algorithm,        'TBLASTX' );
like( $result->algorithm_version, qr/^2\.1\.2/ );
is( $result->query_name, 'HUMBETGLOA' );
is( $result->query_description,
    'Human haplotype C4 beta-globin gene, complete cds.' );
is( $result->query_length,                        3002 );
is( $result->get_statistic('kappa'),              0.135 );
is( $result->get_statistic('lambda'),             0.318 );
is( $result->get_statistic('entropy'),            0.401 );
is( $result->get_statistic('dbletters'),          4662239 );
is( $result->get_statistic('dbentries'),          400 );
is( $result->get_statistic('querylength'),        953 );
is( $result->get_statistic('effectivedblength'),  1535279 );
is( $result->get_statistic('effectivespace'),     1463120887 );
is( $result->get_statistic('effectivespaceused'), 1463120887 );
is( $result->get_statistic('T'),                  13 );
is( $result->get_statistic('X1'),                 16 );
is( $result->get_statistic('X1_bits'),            7.3 );
is( $result->get_statistic('X2'),                 0 );
is( $result->get_statistic('X2_bits'),            '0.0' );
is( $result->get_statistic('S1'),                 41 );
is( $result->get_statistic('S1_bits'),            21.7 );
is( $result->get_statistic('S2'),                 53 );
is( $result->get_statistic('S2_bits'),            27.2 );

is( $result->get_statistic('decayconst'), 0.1 );

is( $result->get_parameter('matrix'), 'BLOSUM62' );

@valid = (
    [ 'gb|AE000479.1|AE000479', 10934, 'AE000479', '0.13', 33.6, 67 ],
    [ 'gb|AE000302.1|AE000302', 10264, 'AE000302', '0.61', 31.3, 62 ],
    [ 'gb|AE000277.1|AE000277', 11653, 'AE000277', '0.84', 30.8, 61 ]
);
$count = 0;

while ( $hit = $result->next_hit ) {
    my $d = shift @valid;
    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->bits,      shift @$d );
    is( $hit->raw_score, shift @$d );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    1057 );
            is( $hsp->query->end,      1134 );
            is( $hsp->query->strand,   1 );
            is( $hsp->strand('query'), $hsp->query->strand );
            is( $hsp->hit->end,        5893 );
            is( $hsp->hit->start,      5816 );
            is( $hsp->hit->strand,     -1 );
            is( $hsp->strand('sbjct'), $hsp->subject->strand );
            is( $hsp->length('total'), 26 );

            float_is( $hsp->evalue, 0.13 );
            is( $hsp->score, 67 );
            is( $hsp->bits,  33.6 );
            is( sprintf( "%.2f", $hsp->percent_identity ),        42.31 );
            is( sprintf( "%.4f", $hsp->frac_identical('query') ), '0.4231' );
            is( sprintf( "%.4f", $hsp->frac_identical('hit') ),   '0.4231' );
            is( $hsp->query->frame(),  0 );
            is( $hsp->hit->frame(),    1 );
            is( $hsp->gaps,            0 );
            is( $hsp->query_string,    'SAYWSIFPPLGCWWSTLGPRGSLSPL' );
            is( $hsp->hit_string,      'AAVWALFPPVGSQWGCLASQWRTSPL' );
            is( $hsp->homology_string, '+A W++FPP+G  W  L  +   SPL' );

            # changed to reflect positional ambiguities, note extra flag
            is(
                join( ' ', $hsp->seq_inds( 'query', 'nomatch', 1 ) ),
                '1063-1065 1090-1095 1099-1104 1108-1113 1117-1125'
            );
            is(
                join( ' ', $hsp->seq_inds( 'hit', 'nomatch', 1 ) ),
                '5825-5833 5837-5842 5846-5851 5855-5860 5885-5887'
            );
            is( $hsp->ambiguous_seq_inds, 'query/subject' );
            is( $hsp->n,                  1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( @valid, 0 );

# test for MarkW bug in blastN

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('a_thaliana.blastn')
);

$result = $searchio->next_result;
like($result->algorithm_reference,qr/Gapped BLAST and PSI-BLAST/);
is( $result->rid, '1012577175-3730-28291' );
is( $result->database_name,
'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS, GSS,or phase 0, 1 or 2 HTGS sequences) '
);
is( $result->database_letters, 4677375331 );
is( $result->database_entries, 1083200 );
is( $result->algorithm,        'BLASTN' );
like( $result->algorithm_version, qr/^2\.2\.1/ );
is( $result->query_name,                          '' );
is( $result->query_length,                        60 );
is( $result->get_parameter('gapopen'),            5 );
is( $result->get_parameter('gapext'),             2 );
is( $result->get_parameter('ktup'),               undef );
is( $result->get_statistic('querylength'),        41 );
is( $result->get_statistic('effectivedblength'),  4656794531 );
is( $result->get_statistic('effectivespace'),     190928575771 );
is( $result->get_statistic('effectivespaceused'), 190928575771 );

is( $result->get_statistic('lambda'),  1.37 );
is( $result->get_statistic('kappa'),   0.711 );
is( $result->get_statistic('entropy'), 1.31 );
is( $result->get_statistic('T'),       0 );
is( $result->get_statistic('A'),       30 );
is( $result->get_statistic('X1'),      '6' );
is( $result->get_statistic('X1_bits'), 11.9 );
is( $result->get_statistic('X2'),      15 );
is( $result->get_statistic('X2_bits'), 29.7 );
is( $result->get_statistic('S1'),      12 );
is( $result->get_statistic('S1_bits'), 24.3 );
is( $result->get_statistic('S2'),      17 );
is( $result->get_statistic('S2_bits'), 34.2 );

is( $result->get_statistic('dbentries'), 1083200 );

@valid = (
    [ 'gb|AY052359.1|', 2826, 'AY052359', '3e-18', 95.6, 48, 1, 60, '1.0000' ],
    [
        'gb|AC002329.2|AC002329', 76170, 'AC002329', '3e-18', 95.6, 48, 1, 60,
        '1.0000'
    ],
    [
        'gb|AF132318.1|AF132318', 5383, 'AF132318', '0.04', 42.1, 21, 35, 55,
        '0.3500'
    ]
);
$count = 0;

while ( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->bits,      shift @$d );
    is( $hit->raw_score, shift @$d );
    is( $hit->start,     shift @$d );
    is( $hit->end,       shift @$d );
    is( sprintf( "%.4f", $hit->frac_aligned_query ), shift @$d );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    1 );
            is( $hsp->query->end,      60 );
            is( $hsp->query->strand,   1 );
            is( $hsp->hit->start,      154 );
            is( $hsp->hit->end,        212 );
            is( $hsp->hit->strand,     1 );
            is( $hsp->length('total'), 60 );
            float_is( $hsp->evalue, 3e-18 );
            is( $hsp->score, 48 );
            is( $hsp->bits,  95.6 );
            is( sprintf( "%.2f", $hsp->percent_identity ),        96.67 );
            is( sprintf( "%.4f", $hsp->frac_identical('query') ), 0.9667 );
            is( sprintf( "%.4f", $hsp->frac_identical('hit') ),   0.9831 );
            is( $hsp->query->frame(), 0 );
            is( $hsp->hit->frame(),   0 );
            is( $hsp->query->seq_id,  undef );
            is( $hsp->hit->seq_id,    'gb|AY052359.1|' );
            is( $hsp->gaps('query'),  0 );
            is( $hsp->gaps('hit'),    1 );
            is( $hsp->gaps,           1 );
            is( $hsp->query_string,
                'aggaatgctgtttaattggaatcgtacaatggagaatttgacggaaatagaatcaacgat'
            );
            is( $hsp->hit_string,
                'aggaatgctgtttaattggaatca-acaatggagaatttgacggaaatagaatcaacgat'
            );
            is( $hsp->homology_string,
                '|||||||||||||||||||||||  |||||||||||||||||||||||||||||||||||'
            );
            my $aln = $hsp->get_aln;
            is( sprintf( "%.2f", $aln->overall_percentage_identity ), 96.67 );
            is( sprintf( "%.2f", $aln->percentage_identity ),         98.31 );
            is( $hsp->n, 1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( @valid, 0 );

#WU-BlastX test

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('dnaEbsub_ecoli.wublastx')
);

$result = $searchio->next_result;
is(
    $result->algorithm_reference, 'Gish, W. (1996-2000) http://blast.wustl.edu
Gish, Warren and David J. States (1993).  Identification of protein coding
regions by database similarity search.  Nat. Genet. 3:266-72.
'
);
is( $result->database_name,    'ecoli.aa' );
is( $result->database_letters, 1358990 );
is( $result->database_entries, 4289 );
is( $result->algorithm,        'BLASTX' );
like( $result->algorithm_version, qr/^2\.0MP\-WashU/ );
is( $result->query_name, 'gi|142864|gb|M10040.1|BACDNAE' );
is( $result->query_description,
    'B.subtilis dnaE gene encoding DNA primase, complete cds' );
is( $result->query_accession,         'M10040.1' );
is( $result->query_gi,                142864 );
is( $result->query_length,            2001 );
is( $result->get_parameter('matrix'), 'blosum62' );

is( $result->get_statistic('lambda'),  0.318 );
is( $result->get_statistic('kappa'),   0.135 );
is( $result->get_statistic('entropy'), 0.401 );

is( $result->get_statistic('dbentries'), 4289 );

@valid = ( [ 'gi|1789447|gb|AAC76102.1|', 581, 'AAC76102', '1.1e-74', 671 ] );
$count = 0;

while ( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->raw_score, shift @$d );
    is( sprintf( "%.4f", $hit->frac_identical('query') ), '0.3640' );
    is( sprintf( "%.4f", $hit->frac_identical('hit') ),   '0.3660' );
    is( sprintf( "%.4f", $hit->frac_conserved('query') ), '0.5370' );
    is( sprintf( "%.4f", $hit->frac_conserved('hit') ),   '0.5400' );
    is( sprintf( "%.4f", $hit->frac_aligned_query ),      '0.6200' );
    is( sprintf( "%.4f", $hit->frac_aligned_hit ),        '0.7100' );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    21 );
            is( $hsp->query->end,      1265 );
            is( $hsp->query->strand,   1 );
            is( $hsp->hit->start,      1 );
            is( $hsp->hit->end,        413 );
            is( $hsp->hit->strand,     0 );
            is( $hsp->length('total'), 421 );
            float_is( $hsp->evalue, 1.1e-74 );
            float_is( $hsp->pvalue, '1.1e-74' );
            is( $hsp->score, 671 );
            is( $hsp->bits,  265.8 );
            is( sprintf( "%.2f", $hsp->percent_identity ), 35.87 );

            is( sprintf( "%.4f", $hsp->frac_identical('query') ), 0.3639 );
            is( sprintf( "%.4f", $hsp->frac_identical('hit') ),   0.3656 );
            is( sprintf( "%.4f", $hsp->frac_conserved('query') ), 0.5373 );
            is( sprintf( "%.2f", $hsp->frac_conserved('hit') ),   0.54 );

            is( sprintf( "%.4f", $hsp->frac_identical('hsp') ), 0.3587 );
            is( sprintf( "%.4f", $hsp->frac_conserved('hsp') ), 0.5297 );

            is( $hsp->query->frame(), 2 );
            is( $hsp->hit->frame(),   0 );
            is( $hsp->gaps('query'),  6 );
            is( $hsp->gaps('hit'),    8 );
            is( $hsp->gaps,           14 );
            is( $hsp->query_string,
'MGNRIPDEIVDQVQKSADIVEVIGDYVQLKKQGRNYFGLCPFHGESTPSFSVSPDKQIFHCFGCGAGGNVFSFLRQMEGYSFAESVSHLADKYQIDFPDDITVHSGARP---ESSGEQKMAEAHELLKKFYHHLLINTKEGQEALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGYFDRFRNRVMFPIHDHHGAVVAFSGRALGSQQPKYMNSPETPLFHKSKLLYNFYKARLHIRKQERAVLFEGFADVYTAVSSDVKESIATMGTSLTDDHVKILRRNVEEIILCYDSDKAGYEATLKASELL---QKKGCKVRVAMIPDGLDPDDYIKKFGGEKFKNDIIDASVTVMAFKMQYFRKGKNLSDEGDRLAYIKDVLKEISTLSGSLEQEVYVKQ'
            );
            is( $hsp->hit_string,
'MAGRIPRVFINDLLARTDIVDLIDARVKLKKQGKNFHACCPFHNEKTPSFTVNGEKQFYHCFGCGAHGNAIDFLMNYDKLEFVETVEELAAMHNLEVPFE----AGSGPSQIERHQRQTLYQLMDGLNTFYQQSL-QQPVATSARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY-DRFRERVMFPIRDKRGRVIGFGGRVLGNDTPKYLNSPETDIFHKGRQLYGLYEAQQDNAEPNRLLVVEGYMDVVALAQYGINYAVASLGTSTTADHIQLLFRATNNVICCYDGDRAGRDAAWRALETALPYMTDGRQLRFMFLPDGEDPDTLVRKEGKEAFEARM-EQAMPLSAFLFNSLMPQVDLSTPDGRARLSTLALPLISQVPGETLR-IYLRQ'
            );
            is( $hsp->homology_string,
'M  RIP   ++ +    DIV++I   V+LKKQG+N+   CPFH E TPSF+V+ +KQ +HCFGCGA GN   FL   +   F E+V  LA  + ++ P +    +G+ P   E    Q + +  + L  FY   L        A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y DRFR RVMFPI D  G V+ F GR LG+  PKY+NSPET +FHK + LY  Y+A+    +  R ++ EG+ DV       +  ++A++GTS T DH+++L R    +I CYD D+AG +A  +A E        G ++R   +PDG DPD  ++K G E F+  + + ++ + AF         +LS    R       L  IS + G   + +Y++Q'
            );
            is(
                join( ' ', $hsp->seq_inds( 'query', 'nomatch', 1 ) ),
'24-29 39-47 54-56 60-71 90-98 129-137 150-152 156-158 180-182 192-194 219-221 228-236 243-251 255-263 267-269 279-284 291-296 300-302 309-311 315-317 321-332 342-344 351-362 366-368 372-374 378-383 387-389 393-398 405-413 417-440 444-449 456-461 468-470 474-476 486-491 495-497 510-518 525-527 531-533 537-557 561-569 573-578 594-599 603-605 609-614 618-620 633-635 654-656 660-665 669-671 678-680 684-686 693-695 705-710 738-740 753-755 759-761 768-773 786-797 801-806 810-812 819-821 831-833 840-860 864-869 894-896 900-902 921-923 927-938 945-947 957-959 972-974 981-986 993-995 999-1013 1017-1019 1029-1037 1050-1052 1062-1067 1077-1079 1083-1085 1089-1091 1098-1103 1107-1109 1113-1115 1122-1124 1128-1130 1137-1163 1173-1184 1188-1208 1212-1217 1224-1226 1230-1232 1236-1244 1248-1250'
            );
            is(
                join( ' ', $hsp->seq_inds( 'query', 'mismatch', 1 ) ),
'24-29 39-47 54-56 60-71 90-98 129-137 150-152 156-158 180-182 192-194 219-221 228-236 243-251 255-263 267-269 279-284 291-296 300-302 309-311 315-317 342-344 351-362 366-368 372-374 378-383 387-389 393-398 405-413 420-440 444-449 456-461 468-470 474-476 486-491 495-497 510-518 525-527 531-533 537-557 561-569 573-578 594-599 603-605 609-614 633-635 654-656 660-665 669-671 678-680 684-686 693-695 705-710 738-740 753-755 759-761 768-773 786-797 801-806 810-812 819-821 831-833 840-860 864-869 894-896 900-902 921-923 927-938 945-947 957-959 972-974 981-986 993-995 999-1013 1017-1019 1029-1037 1050-1052 1062-1067 1077-1079 1083-1085 1089-1091 1098-1103 1113-1115 1122-1124 1128-1130 1137-1163 1173-1184 1188-1208 1212-1217 1224-1226 1230-1232 1236-1244'
            );
            is(
                join( ' ', $hsp->seq_inds( 'hit', 'nomatch', 1 ) ),
'2 3 7-9 12 14-17 24-26 37-39 44 46 54 58 67 70-72 75-77 79-81 83 87 88 91 92 94 97 99 104 106-108 110-113 115 117 119 120 122 124 125 128-130 132-138 140 141 144 145 148 150 154 155 157 162-164 167 169 171-177 179-181 183 184 190 191 193 195 196 202 209 211 212 214 217 219 222 226 227 237 242 244 247 248 253-256 258 259 261 264 268 271-277 279 280 289 291 298 300-303 306 310 315 318 319 322 324-331 333 337-339 344 348 349 353 355 357 360 361 364 367 369 372-380 384-387 389-395 397 398 401 403 405-407'
            );
            is(
                join( ' ', $hsp->seq_inds( 'hit', 'mismatch', 1 ) ),
'2 3 7-9 12 14-17 24-26 37-39 44 46 54 58 67 70-72 75-77 79-81 83 87 88 91 92 94 97 99 104 110-113 115 117 119 120 122 124 125 128-130 132-138 140 141 144 145 148 150 154 155 157 162-164 167 169 171-177 179-181 183 184 190 191 193 195 196 202 209 211 212 214 217 219 222 226 227 237 242 244 247 248 253-256 258 259 261 264 268 271-277 279 280 289 291 298 300-303 306 310 315 318 319 322 324 325 329-331 333 337-339 344 348 349 353 355 357 360 361 364 367 369 372-380 384-387 389-395 397 398 401 403 405-407'
            );
            is( join( ' ', $hsp->seq_inds( 'query', 'gaps', 1 ) ), '347 1004' );
            is( join( ' ', $hsp->seq_inds( 'hit', 'gaps', 1 ) ),
                '100 131 197 362 408' );
            is( $hsp->ambiguous_seq_inds, 'query' );
            is( $hsp->n,                  1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( @valid, 0 );

#Trickier WU-Blast
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('tricky.wublast')
);
$result = $searchio->next_result;
my $hits_left = 1;
while ( my $hit = $result->next_hit ) {

# frac_aligned_hit used to be over 1, frac_identical & frac_conserved are still too wrong
  TODO: {
        local $TODO = 'frac_identical & frac_conserved are still too wrong';
        cmp_ok sprintf( "%.3f", $hit->frac_identical ), '>',  0.9;
        cmp_ok sprintf( "%.3f", $hit->frac_conserved ), '<=', 1;
    }
    is( sprintf( "%.2f", $hit->frac_aligned_query ), '0.92' );
    is( sprintf( "%.2f", $hit->frac_aligned_hit ),   '0.91' );
    $hits_left--;
}
is( $hits_left, 0 );

# More frac_ method testing, this time on ncbi blastn
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('frac_problems.blast')
);
my @expected = ( "1.000", "0.943" );
while ( my $result = $searchio->next_result ) {
    my $hit = $result->next_hit;
    is( $hit->frac_identical, shift @expected );
}
is( @expected, 0 );

# And even more: frac_aligned_query should never be over 1!
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('frac_problems2.blast')
);
$result = $searchio->next_result;
$hit    = $result->next_hit;
is $hit->frac_aligned_query, 0.97;

# Also, start and end should be sane
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('frac_problems3.blast')
);
$result = $searchio->next_result;
$hit    = $result->next_hit;
is $hit->start('sbjct'), 207;
is $hit->end('sbjct'),   1051;

#WU-TBlastN test

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('dnaEbsub_ecoli.wutblastn')
);

$result = $searchio->next_result;
is(
    $result->algorithm_reference, 'Gish, W. (1996-2000) http://blast.wustl.edu
'
);
is( $result->database_name,    'ecoli.nt' );
is( $result->database_letters, 4662239 );
is( $result->database_entries, 400 );
is( $result->algorithm,        'TBLASTN' );
like( $result->algorithm_version, qr/^2\.0MP\-WashU/ );
is( $result->query_name,              'gi|142865|gb|AAA22406.1|' );
is( $result->query_description,       'DNA primase' );
is( $result->query_accession,         'AAA22406.1' );
is( $result->query_gi,                142865 );
is( $result->query_length,            603 );
is( $result->get_parameter('matrix'), 'blosum62' );

is( $result->get_statistic('lambda'),  '0.320' );
is( $result->get_statistic('kappa'),   0.136 );
is( $result->get_statistic('entropy'), 0.387 );

is( $result->get_statistic('dbentries'), 400 );

@valid =
  ( [ 'gi|1789441|gb|AE000388.1|AE000388', 10334, 'AE000388', '1.4e-73', 671 ]
  );
$count = 0;

while ( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->raw_score, shift @$d );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    1 );
            is( $hsp->query->end,      415 );
            is( $hsp->query->strand,   0 );
            is( $hsp->hit->start,      4778 );
            is( $hsp->hit->end,        6016 );
            is( $hsp->hit->strand,     1 );
            is( $hsp->length('total'), 421 );
            float_is( $hsp->evalue, 1.4e-73 );
            float_is( $hsp->pvalue, 1.4e-73 );
            is( $hsp->score, 671 );
            is( $hsp->bits,  265.8 );
            is( sprintf( "%.2f", $hsp->percent_identity ),        35.87 );
            is( sprintf( "%.4f", $hsp->frac_identical('hit') ),   0.3656 );
            is( sprintf( "%.4f", $hsp->frac_identical('query') ), 0.3639 );
            is( sprintf( "%.4f", $hsp->frac_conserved('hsp') ),   0.5297 );
            is( $hsp->query->frame(), 0 );
            is( $hsp->hit->frame(),   1 );
            is( $hsp->gaps('query'),  6 );
            is( $hsp->gaps('hit'),    8 );
            is( $hsp->gaps,           14 );
            is( $hsp->query_string,
'MGNRIPDEIVDQVQKSADIVEVIGDYVQLKKQGRNYFGLCPFHGESTPSFSVSPDKQIFHCFGCGAGGNVFSFLRQMEGYSFAESVSHLADKYQIDFPDDITVHSGARP---ESSGEQKMAEAHELLKKFYHHLLINTKEGQEALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGYFDRFRNRVMFPIHDHHGAVVAFSGRALGSQQPKYMNSPETPLFHKSKLLYNFYKARLHIRKQERAVLFEGFADVYTAVSSDVKESIATMGTSLTDDHVKILRRNVEEIILCYDSDKAGYEATLKASELL---QKKGCKVRVAMIPDGLDPDDYIKKFGGEKFKNDIIDASVTVMAFKMQYFRKGKNLSDEGDRLAYIKDVLKEISTLSGSLEQEVYVKQ'
            );
            is( $hsp->hit_string,
'MAGRIPRVFINDLLARTDIVDLIDARVKLKKQGKNFHACCPFHNEKTPSFTVNGEKQFYHCFGCGAHGNAIDFLMNYDKLEFVETVEELAAMHNLEVPFE----AGSGPSQIERHQRQTLYQLMDGLNTFYQQSL-QQPVATSARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY-DRFRERVMFPIRDKRGRVIGFGGRVLGNDTPKYLNSPETDIFHKGRQLYGLYEAQQDNAEPNRLLVVEGYMDVVALAQYGINYAVASLGTSTTADHIQLLFRATNNVICCYDGDRAGRDAAWRALETALPYMTDGRQLRFMFLPDGEDPDTLVRKEGKEAFEARM-EQAMPLSAFLFNSLMPQVDLSTPDGRARLSTLALPLISQVPGETLR-IYLRQ'
            );
            is( $hsp->homology_string,
'M  RIP   ++ +    DIV++I   V+LKKQG+N+   CPFH E TPSF+V+ +KQ +HCFGCGA GN   FL   +   F E+V  LA  + ++ P +    +G+ P   E    Q + +  + L  FY   L        A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y DRFR RVMFPI D  G V+ F GR LG+  PKY+NSPET +FHK + LY  Y+A+    +  R ++ EG+ DV       +  ++A++GTS T DH+++L R    +I CYD D+AG +A  +A E        G ++R   +PDG DPD  ++K G E F+  + + ++ + AF         +LS    R       L  IS + G   + +Y++Q'
            );
            is(
                join( ' ', $hsp->seq_inds( 'query', 'nomatch', 1 ) ),
'2 3 7-9 12 14-17 24-26 37-39 44 46 54 58 67 70-72 75-77 79-81 83 87 88 91 92 94 97 99 101-104 108 111-114 116 118 120 121 123 125 126 129-131 133-140 142 143 146 147 150 152 156 157 159 164-166 169 171 173-179 181-183 185 186 192 193 195 197 198 200 205 212 214 215 217 220 222 225 229 230 240 245 247 250 251 256-259 261 262 264 267 271 274-280 282 283 292 294 301 303-306 309 313 318 321 322 325 327-331 333 337-339 344 348 349 353 355 357 360 361 363 365 368 370 373-381 385-388 390-396 398 399 402 404 406-408 410'
            );
            is(
                join( ' ', $hsp->seq_inds( 'hit', 'nomatch', 1 ) ),
'4781-4786 4796-4804 4811-4813 4817-4828 4847-4855 4886-4894 4907-4909 4913-4915 4937-4939 4949-4951 4976-4978 4985-4993 5000-5008 5012-5020 5024-5026 5036-5041 5048-5053 5057-5059 5066-5068 5072-5074 5087-5089 5093-5101 5105-5116 5120-5122 5126-5128 5132-5137 5141-5143 5147-5152 5159-5167 5171-5191 5195-5200 5207-5212 5219-5221 5225-5227 5237-5242 5246-5248 5261-5269 5276-5278 5282-5284 5288-5308 5312-5320 5324-5329 5345-5350 5354-5356 5360-5365 5381-5383 5402-5404 5408-5413 5417-5419 5426-5428 5432-5434 5441-5443 5453-5458 5486-5488 5501-5503 5507-5509 5516-5521 5534-5545 5549-5554 5558-5560 5567-5569 5579-5581 5588-5608 5612-5617 5642-5644 5648-5650 5669-5671 5675-5686 5693-5695 5705-5707 5720-5722 5729-5734 5741-5743 5747-5770 5774-5776 5786-5794 5807-5809 5819-5824 5834-5836 5840-5842 5846-5848 5855-5860 5867-5869 5876-5878 5882-5884 5891-5917 5927-5938 5942-5962 5966-5971 5978-5980 5984-5986 5990-5998'
            );
            is( join( ' ', $hsp->seq_inds( 'query', 'gaps', 1 ) ), '109 328' );
            is( join( ' ', $hsp->seq_inds( 'hit', 'gaps', 1 ) ),
                '5077 5170 5368 5863 6001' );
            is( $hsp->ambiguous_seq_inds, 'subject' );
            is( $hsp->n,                  1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( $count, 1 );

# WU-BLAST TBLASTX
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('dnaEbsub_ecoli.wutblastx')
);

$result = $searchio->next_result;
is(
    $result->algorithm_reference, 'Gish, W. (1996-2000) http://blast.wustl.edu
'
);
is( $result->database_name,    'ecoli.nt' );
is( $result->database_letters, 4662239 );
is( $result->database_entries, 400 );
is( $result->algorithm,        'TBLASTX' );
like( $result->algorithm_version, qr/^2\.0MP\-WashU/ );
is( $result->query_name, 'gi|142864|gb|M10040.1|BACDNAE' );
is( $result->query_description,
    'B.subtilis dnaE gene encoding DNA primase, complete cds' );
is( $result->query_accession,         'M10040.1' );
is( $result->query_gi,                142864 );
is( $result->query_length,            2001 );
is( $result->get_parameter('matrix'), 'blosum62' );

is( $result->get_statistic('lambda'),    0.318 );
is( $result->get_statistic('kappa'),     0.135 );
is( $result->get_statistic('entropy'),   0.401 );
is( $result->get_statistic('dbentries'), 400 );

@valid = (
    [
        'gi|1789441|gb|AE000388.1|AE000388',
        10334, 'AE000388', '6.4e-70', 318, 148.6
    ],
    [ 'gi|2367383|gb|AE000509.1|AE000509', 10589, 'AE000509', 1, 59, 29.9 ]
);
$count = 0;

while ( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );

    # using e here to deal with 0.9992 coming out right here as well
    float_is( $hit->significance, shift @$d );
    is( $hit->raw_score, shift @$d );
    is( $hit->bits,      shift @$d );
    if ( $count == 0 ) {
        my $hspcounter = 0;
        while ( my $hsp = $hit->next_hsp ) {
            $hspcounter++;
            if ( $hspcounter == 3 ) {

                # let's actually look at the 3rd HSP
                is( $hsp->query->start,    441 );
                is( $hsp->query->end,      617 );
                is( $hsp->query->strand,   1 );
                is( $hsp->hit->start,      5192 );
                is( $hsp->hit->end,        5368 );
                is( $hsp->hit->strand,     1 );
                is( $hsp->length('total'), 59 );
                float_is( $hsp->evalue, 6.4e-70 );
                float_is( $hsp->pvalue, 6.4e-70 );
                is( $hsp->score, 85 );
                is( $hsp->bits,  41.8 );
                is( sprintf( "%.2f", $hsp->percent_identity ),        '32.20' );
                is( sprintf( "%.3f", $hsp->frac_identical('hit') ),   0.322 );
                is( sprintf( "%.3f", $hsp->frac_identical('query') ), 0.322 );
                is( sprintf( "%.4f", $hsp->frac_conserved('hsp') ),   0.4746 );
                is( $hsp->query->frame(), 2 );
                is( $hsp->hit->frame(),   1 );
                is( $hsp->gaps('query'),  0 );
                is( $hsp->gaps('hit'),    0 );
                is( $hsp->gaps,           0 );
                is( $hsp->n,              1 );
                is( $hsp->query_string,
'ALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGY'
                );
                is( $hsp->hit_string,
'ARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY'
                );
                is( $hsp->homology_string,
'A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y'
                );
                is(
                    join( ' ', $hsp->seq_inds( 'query', 'nomatch', 1 ) ),
'444-449 456-461 468-470 474-476 486-491 495-497 510-518 525-527 531-533 537-557 561-569 573-578 594-599 603-605 609-614'
                );
                is(
                    join( ' ', $hsp->seq_inds( 'hit', 'nomatch', 1 ) ),
'5195-5200 5207-5212 5219-5221 5225-5227 5237-5242 5246-5248 5261-5269 5276-5278 5282-5284 5288-5308 5312-5320 5324-5329 5345-5350 5354-5356 5360-5365'
                );
                is( $hsp->ambiguous_seq_inds, 'query/subject' );
                last;
            }
        }
        is( $hspcounter, 3 );
    }
    elsif ( $count == 1 ) {
        my $hsps_to_do = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    587 );
            is( $hsp->query->end,      706 );
            is( $hsp->query->strand,   -1 );
            is( $hsp->hit->start,      4108 );
            is( $hsp->hit->end,        4227 );
            is( $hsp->hit->strand,     -1 );
            is( $hsp->length('total'), 40 );
            float_is( $hsp->evalue, 7.1 );
            float_is( $hsp->pvalue, '1.00' );
            is( $hsp->score, 59 );
            is( $hsp->bits,  29.9 );
            is( sprintf( "%.2f", $hsp->percent_identity ),        '37.50' );
            is( sprintf( "%.4f", $hsp->frac_identical('hit') ),   '0.3750' );
            is( sprintf( "%.4f", $hsp->frac_identical('query') ), '0.3750' );
            is( sprintf( "%.4f", $hsp->frac_conserved('hsp') ),   '0.4750' );
            is( $hsp->query->frame(), 2 );
            is( $hsp->hit->frame(),   2 );
            is( $hsp->gaps('query'),  0 );
            is( $hsp->gaps('hit'),    0 );
            is( $hsp->gaps,           0 );
            is( $hsp->n,              1 );
            is( $hsp->query_string,
                'WLPRALPEKATTAP**SWIGNMTRFLKRSKYPLPSSRLIR' );
            is( $hsp->hit_string, 'WLSRTTVGSSTVSPRTFWITRMKVKLSSSKVTLPSTKSTR' );
            is( $hsp->homology_string,
                'WL R     +T +P   WI  M   L  SK  LPS++  R' );
            $hsps_to_do--;
            last;
        }
        is( $hsps_to_do, 0 );
    }
    last if ( $count++ > @valid );
}
is( $count, 2 );

# WU-BLAST -echofilter option test (Bug 2388)
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('echofilter.wublastn')
);

$result = $searchio->next_result;
is(
    $result->algorithm_reference, 'Gish, W. (1996-2006) http://blast.wustl.edu
'
);
is( $result->database_name,    'NM_003201.fa' );
is( $result->database_letters, 1936 );
is( $result->database_entries, 1 );
is( $result->algorithm,        'BLASTN' );
like( $result->algorithm_version, qr/^2\.0MP\-WashU/ );
like( $result->query_name,
qr/ref|NM_003201.1| Homo sapiens transcription factor A, mitochondrial \(TFAM\), mRNA/
);
is( $result->query_accession, 'NM_003201.1' );

is( $result->query_length,               1936 );
is( $result->get_statistic('lambda'),    0.192 );
is( $result->get_statistic('kappa'),     0.182 );
is( $result->get_statistic('entropy'),   0.357 );
is( $result->get_statistic('dbletters'), 1936 );
is( $result->get_statistic('dbentries'), 1 );
is( $result->get_parameter('matrix'),    '+5,-4' );

@valid = ( [ 'ref|NM_003201.1|', 1936, 'NM_003201', '0', 9680 ], );
$count = 0;
while ( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->raw_score, shift @$d );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    1 );
            is( $hsp->query->end,      1936 );
            is( $hsp->hit->start,      1 );
            is( $hsp->hit->end,        1936 );
            is( $hsp->length('total'), 1936 );

            float_is( $hsp->evalue, 0. );
            float_is( $hsp->pvalue, '0.' );
            is( $hsp->score,                   9680 );
            is( $hsp->bits,                    1458.4 );
            is( $hsp->percent_identity,        100 );
            is( $hsp->frac_identical('query'), 1.00 );
            is( $hsp->frac_identical('hit'),   1.00 );
            is( $hsp->gaps,                    0 );
            is( $hsp->n,                       1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( @valid, 0 );

# Do a multiblast report test
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('multi_blast.bls')
);

@expected = qw(CATH_RAT CATL_HUMAN CATL_RAT PAPA_CARPA);
my $results_left = 4;
while ( my $result = $searchio->next_result ) {
    like($result->algorithm_reference, qr/Gapped BLAST and PSI-BLAST/);
    is( $result->query_name, shift @expected, "Multiblast query test" );
    $results_left--;
}
is( $results_left, 0 );

# Test GCGBlast parsing

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('test.gcgblast')
);
$result = $searchio->next_result();
like($result->algorithm_reference,qr/Gapped BLAST and PSI-BLAST/);
is( $result->query_name,                   '/v0/people/staji002/test.gcg' );
is( $result->algorithm,                    'BLASTP' );
is( $result->algorithm_version,            '2.2.1 [Apr-13-2001]' );
is( $result->database_name,                'pir' );
is( $result->database_entries,             274514 );
is( $result->database_letters,             93460074 );
is( $result->get_statistic('querylength'), 44 );
is( $result->get_statistic('effectivedblength'),  65459646 );
is( $result->get_statistic('effectivespace'),     2880224424 );
is( $result->get_statistic('effectivespaceused'), 2880224424 );

$hit = $result->next_hit;
is( $hit->description, 'F22B7.10 protein - Caenorhabditis elegans' );
is( $hit->name,        'PIR2:S44629' );
is( $hit->length,      628 );
is( $hit->accession,   'PIR2:S44629' );
float_is( $hit->significance, 2e-08 );
is( $hit->raw_score, 136 );
is( $hit->bits,      '57.0' );
$hsp = $hit->next_hsp;
float_is( $hsp->evalue, 2e-08 );
is( $hsp->bits,                    '57.0' );
is( $hsp->score,                   136 );
is( int( $hsp->percent_identity ), 28 );
is( sprintf( "%.2f", $hsp->frac_identical('query') ), 0.29 );
is( $hsp->frac_conserved('total'), 69 / 135 );
is( $hsp->gaps('total'),           8 );
is( $hsp->gaps('hit'),             6 );
is( $hsp->gaps('query'),           2 );

is( $hsp->hit->start,   342 );
is( $hsp->hit->end,     470 );
is( $hsp->query->start, 3 );
is( $hsp->query->end,   135 );

is( $hsp->query_string,
'CAAEFDFMEKETPLRYTKTXXXXXXXXXXXXXXRKIISDMWGVLAKQQTHVRKHQFDHGELVYHALQLLAYTALGILIMRLKLFLTPYMCVMASLICSRQLFGW--LFCKVHPGAIVFVILAAMSIQGSANLQTQ'
);
is( $hsp->hit_string,
'CSAEFDFIQYSTIEKLCGTLLIPLALISLVTFVFNFVKNT-NLLWRNSEEIG----ENGEILYNVVQLCCSTVMAFLIMRLKLFMTPHLCIVAALFANSKLLGGDRISKTIRVSALVGVI-AILFYRGIPNIRQQ'
);
is( $hsp->homology_string,
'C+AEFDF++  T  +   T                 + +   +L +    +     ++GE++Y+ +QL   T +  LIMRLKLF+TP++C++A+L  + +L G   +   +   A+V VI A +  +G  N++ Q'
);

#test all the database accession number formats
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('testdbaccnums.out')
);
$result = $searchio->next_result;
like($result->algorithm_reference,qr/Gapped BLAST and PSI-BLAST/);
is( $result->rid,                                 '1036160600-011802-21377' );
is( $result->get_statistic('querylength'),        9 );
is( $result->get_statistic('effectivedblength'),  35444647 );
is( $result->get_statistic('effectivespace'),     319001823 );
is( $result->get_statistic('effectivespaceused'), 319001823 );

@valid = (
    [ 'pir||T14789',           'T14789',    'T14789', 'CAB53709', 'AAH01726' ],
    [ 'gb|NP_065733.1|CYT19',  'NP_065733', 'CYT19' ],
    [ 'emb|XP_053690.4|Cyt19', 'XP_053690' ],
    [ 'dbj|NP_056277.2|DKFZP586L0724',      'NP_056277' ],
    [ 'prf||XP_064862.2',                   'XP_064862' ],
    [ 'pdb|BAB13968.1|1',                   'BAB13968' ],
    [ 'sp|Q16478|GLK5_HUMAN',               'Q16478' ],
    [ 'pat|US|NP_002079.2',                 'NP_002079' ],
    [ 'bbs|NP_079463.2|',                   'NP_079463' ],
    [ 'gnl|db1|NP_002444.1',                'NP_002444' ],
    [ 'ref|XP_051877.1|',                   'XP_051877' ],
    [ 'lcl|AAH16829.1|',                    'AAH16829' ],
    [ 'gi|1|gb|NP_065733.1|CYT19',          'NP_065733' ],
    [ 'gi|2|emb|XP_053690.4|Cyt19',         'XP_053690' ],
    [ 'gi|3|dbj|NP_056277.2|DKFZP586L0724', 'NP_056277' ],
    [ 'gi|4|pir||T14789',                   'T14789' ],
    [ 'gi|5|prf||XP_064862.2',              'XP_064862' ],
    [ 'gi|6|pdb|BAB13968.1|1',              'BAB13968' ],
    [ 'gi|7|sp|Q16478|GLK5_HUMAN',          'Q16478' ],
    [ 'gi|8|pat|US|NP_002079.2',            'NP_002079' ],
    [ 'gi|9|bbs|NP_079463.2|',              'NP_079463' ],
    [ 'gi|10|gnl|db1|NP_002444.1',          'NP_002444' ],
    [ 'gi|11|ref|XP_051877.1|',             'XP_051877' ],
    [ 'gi|12|lcl|AAH16829.1|',              'AAH16829' ],
    [ 'MY_test_ID',                         'MY_test_ID' ]
);

$hit = $result->next_hit;
my $d = shift @valid;
is( $hit->name,      shift @$d );
is( $hit->accession, shift @$d );
my @accnums = $hit->each_accession_number;
foreach my $a (@accnums) {
    is( $a, shift @$d );
}
$d   = shift @valid;
$hit = $result->next_hit;
is( $hit->name,      shift @$d );
is( $hit->accession, shift @$d );
is( $hit->locus,     shift @$d );

$hits_left = 23;
while ( $hit = $result->next_hit ) {
    my $d = shift @valid;
    is( $hit->name,      shift @$d );
    is( $hit->accession, shift @$d );
    $hits_left--;
}
is( $hits_left, 0 );

# Parse MEGABLAST

# parse the BLAST-like output
my $infile = test_input_file('503384.MEGABLAST.2');
my $in     = Bio::SearchIO->new(
    -file   => $infile,
    -format => 'blast'
);    # this is megablast blast-like output
my $r        = $in->next_result;
my @dcompare = (
    [
        'Contig3700', 5631, 396, 785, '0.0', 785, '0.0', 396, 639, 12, 8723,
        9434, 1, 4083, 4794, -1
    ],
    [
        'Contig3997', 12734, 335, 664, '0.0', 664, '0.0', 335, 401, 0, 1282,
        1704, 1, 1546, 1968, -1
    ],
    [
        'Contig634', 858, 245, 486, '1e-136', 486,
        '1e-136',    245, 304, 3,   7620,     7941,
        1,           1,   321, -1
    ],
    [
        'Contig1853', 2314, 171,  339, '1e-91', 339,
        '1e-91',      171,  204,  0,   6406,    6620,
        1,            1691, 1905, 1
    ]
);

like($r->algorithm_reference,qr/A greedy algorithm for aligning DNA sequences/);
is( $r->algorithm,                           'MEGABLAST' );
is( $r->query_name,                          '503384' );
is( $r->query_description,                   '11337 bp 2 contigs' );
is( $r->query_length,                        11337 );
is( $r->database_name,                       'cneoA.nt' );
is( $r->database_letters,                    17206226 );
is( $r->database_entries,                    4935 );
is( $r->get_statistic('querylength'),        11318 );
is( $r->get_statistic('effectivedblength'),  17112461 );
is( $r->get_statistic('effectivespace'),     193678833598 );
is( $r->get_statistic('effectivespaceused'), 0 );

$hits_left = 4;
while ( my $hit = $r->next_hit ) {
    my $d = shift @dcompare;
    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->raw_score, shift @$d );
    is( $hit->bits,      shift @$d );
    float_is( $hit->significance, shift @$d );

    my $hsp = $hit->next_hsp;
    is( $hsp->bits, shift @$d );
    float_is( $hsp->evalue, shift @$d );
    is( $hsp->score,         shift @$d );
    is( $hsp->num_identical, shift @$d );
    is( $hsp->gaps('total'), shift @$d );
    is( $hsp->query->start,  shift @$d );
    is( $hsp->query->end,    shift @$d );
    is( $hsp->query->strand, shift @$d );
    is( $hsp->hit->start,    shift @$d );
    is( $hsp->hit->end,      shift @$d );
    is( $hsp->hit->strand,   shift @$d );
    is( $hsp->n,             1 );
    $hits_left--;
}
is( $hits_left, 0 );

# Let's test RPS-BLAST

my $parser = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('ecoli_domains.rpsblast')
);

$r = $parser->next_result;
is( $r->algorithm,                           'RPS-BLAST(BLASTP)');
is( $r->algorithm_version,                   '2.2.4 [Aug-26-2002]');
is( $r->algorithm_reference,                 undef );
is( $r->query_name,                          'gi|1786183|gb|AAC73113.1|' );
is( $r->query_gi,                            1786183 );
is( $r->num_hits,                            7 );
is( $r->get_statistic('querylength'),        438 );
is( $r->get_statistic('effectivedblength'),  31988 );
is( $r->get_statistic('effectivespace'),     14010744 );
is( $r->get_statistic('effectivespaceused'), 24054976 );
$hit = $r->next_hit;
is( $hit->name, 'gnl|CDD|3919' );
float_is( $hit->significance, 0.064 );
is( $hit->bits,      28.3 );
is( $hit->raw_score, 63 );
$hsp = $hit->next_hsp;
is( $hsp->query->start, 599 );
is( $hsp->query->end,   655 );
is( $hsp->hit->start,   23 );
is( $hsp->hit->end,     76 );

# Test PSI-BLAST parsing

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('psiblastreport.out')
);

$result = $searchio->next_result;
like($result->algorithm_reference, qr/Gapped BLAST and PSI-BLAST/);
is( $result->database_name,    '/home/peter/blast/data/swissprot.pr' );
is( $result->database_entries, 88780 );
is( $result->database_letters, 31984247 );

is( $result->algorithm, 'BLASTP' );
like( $result->algorithm_version, qr/^2\.0\.14/ );
is( $result->query_name,             'CYS1_DICDI' );
is( $result->query_length,           343 );
is( $result->get_statistic('kappa'), 0.0491 );
cmp_ok( $result->get_statistic('lambda'),  '==', 0.270 );
cmp_ok( $result->get_statistic('entropy'), '==', 0.230 );
is( $result->get_statistic('dbletters'),           31984247 );
is( $result->get_statistic('dbentries'),           88780 );
is( $result->get_statistic('effective_hsplength'), 49 );
is( $result->get_statistic('querylength'),         294 );
is( $result->get_statistic('effectivedblength'),   27634027 );
is( $result->get_statistic('effectivespace'),      8124403938 );
is( $result->get_statistic('effectivespaceused'),  8124403938 );
is( $result->get_parameter('matrix'),              'BLOSUM62' );
is( $result->get_parameter('gapopen'),             11 );
is( $result->get_parameter('gapext'),              1 );

my @valid_hit_data = (
    [ 'sp|P04988|CYS1_DICDI', 343, 'P04988', '0',     721 ],
    [ 'sp|P43295|A494_ARATH', 313, 'P43295', '1e-75', 281 ],
    [ 'sp|P25804|CYSP_PEA',   363, 'P25804', '1e-74', 278 ]
);
my @valid_iter_data = (
    [ 127, 127, 0,   109, 18, 0, 0,   0, 0 ],
    [ 157, 40,  117, 2,   38, 0, 109, 3, 5 ]
);
my $iter_count = 0;

while ( $iter = $result->next_iteration ) {
    $iter_count++;
    my $di = shift @valid_iter_data;
    is( $iter->number, $iter_count );

    is( $iter->num_hits,                                shift @$di );
    is( $iter->num_hits_new,                            shift @$di );
    is( $iter->num_hits_old,                            shift @$di );
    is( scalar( $iter->newhits_below_threshold ),       shift @$di );
    is( scalar( $iter->newhits_not_below_threshold ),   shift @$di );
    is( scalar( $iter->newhits_unclassified ),          shift @$di );
    is( scalar( $iter->oldhits_below_threshold ),       shift @$di );
    is( scalar( $iter->oldhits_newly_below_threshold ), shift @$di );
    is( scalar( $iter->oldhits_not_below_threshold ),   shift @$di );

    my $hit_count = 0;
    if ( $iter_count == 1 ) {
        while ( $hit = $result->next_hit ) {
            my $d = shift @valid_hit_data;

            is( $hit->name,      shift @$d );
            is( $hit->length,    shift @$d );
            is( $hit->accession, shift @$d );
            float_is( $hit->significance, shift @$d );
            is( $hit->bits, shift @$d );

            if ( $hit_count == 1 ) {
                my $hsps_left = 1;
                while ( my $hsp = $hit->next_hsp ) {
                    is( $hsp->query->start,    32 );
                    is( $hsp->query->end,      340 );
                    is( $hsp->hit->start,      3 );
                    is( $hsp->hit->end,        307 );
                    is( $hsp->length('total'), 316 );
                    is( $hsp->start('hit'),    $hsp->hit->start );
                    is( $hsp->end('query'),    $hsp->query->end );
                    is( $hsp->strand('sbjct'), $hsp->subject->strand )
                      ;    # alias for hit
                    float_is( $hsp->evalue, 1e-75 );
                    is( $hsp->score, 712 );
                    is( $hsp->bits,  281 );
                    is( sprintf( "%.1f", $hsp->percent_identity ), 46.5 );
                    is( sprintf( "%.4f", $hsp->frac_identical('query') ),
                        0.4757 );
                    is( sprintf( "%.3f", $hsp->frac_identical('hit') ), 0.482 );
                    is( $hsp->gaps, 18 );
                    is( $hsp->n,    1 );
                    $hsps_left--;
                }
                is( $hsps_left, 0 );
            }
            last if ( $hit_count++ > @valid_hit_data );
        }
    }
}
is( @valid_hit_data,  0 );
is( @valid_iter_data, 0 );

# Test filtering

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('ecolitst.bls'),
    '-signif' => 1e-100
);

@valid = qw(gb|AAC73113.1|);
$r     = $searchio->next_result;

while ( my $hit = $r->next_hit ) {
    is( $hit->name, shift @valid );
}

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('ecolitst.bls'),
    '-score'  => 100
);

@valid = qw(gb|AAC73113.1| gb|AAC76922.1| gb|AAC76994.1|);
$r     = $searchio->next_result;

while ( my $hit = $r->next_hit ) {
    is( $hit->name, shift @valid );
}
is( @valid, 0 );

$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('ecolitst.bls'),
    '-bits'   => 200
);

@valid = qw(gb|AAC73113.1| gb|AAC76922.1|);
$r     = $searchio->next_result;

while ( my $hit = $r->next_hit ) {
    is( $hit->name, shift @valid );
}
is( @valid, 0 );

my $filt_func = sub {
    my $hit = shift;
    $hit->frac_identical('query') >= 0.31;
};

$searchio = Bio::SearchIO->new(
    '-format'     => 'blast',
    '-file'       => test_input_file('ecolitst.bls'),
    '-hit_filter' => $filt_func
);

@valid = qw(gb|AAC73113.1| gb|AAC76994.1|);
$r     = $searchio->next_result;

while ( my $hit = $r->next_hit ) {
    is( $hit->name, shift @valid );
}
is( @valid, 0 );

# bl2seq parsing testing

# this is blastp bl2seq
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('bl2seq.out')
);
$result = $searchio->next_result;
isa_ok( $result, 'Bio::Search::Result::ResultI' );
is( $result->query_name,                          '' );
is( $result->algorithm,                           'BLASTP' );
is( $result->algorithm_reference,                 undef );
is( $result->get_statistic('querylength'),        320 );
is( $result->get_statistic('effectivedblength'),  339 );
is( $result->get_statistic('effectivespace'),     108480 );
is( $result->get_statistic('effectivespaceused'), 108480 );
$hit = $result->next_hit;
is( $hit->name,   'ALEU_HORVU' );
is( $hit->length, 362 );
$hsp = $hit->next_hsp;
is( $hsp->score,                481 );
is( $hsp->bits,                 191 );
is( int $hsp->percent_identity, 34 );
float_is( $hsp->evalue, 2e-53 );
is( int( $hsp->frac_conserved * $hsp->length ), 167 );
is( $hsp->query->start,                         28 );
is( $hsp->query->end,                           343 );
is( $hsp->hit->start,                           60 );
is( $hsp->hit->end,                             360 );
is( $hsp->gaps,                                 27 );

# this is blastn bl2seq
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('bl2seq.blastn.rev')
);
$result = $searchio->next_result;
isa_ok( $result, 'Bio::Search::Result::ResultI' );
is( $result->query_name,                          '' );
is( $result->algorithm,                           'BLASTN' );
is( $result->algorithm_reference,                 undef );
is( $result->query_length,                        180 );
is( $result->get_statistic('querylength'),        174 );
is( $result->get_statistic('effectivedblength'),  173 );
is( $result->get_statistic('effectivespace'),     30102 );
is( $result->get_statistic('effectivespaceused'), 30102 );
$hit = $result->next_hit;
is( $hit->length, 179 );
is( $hit->name,   'human' );
$hsp = $hit->next_hsp;
is( $hsp->score,                27 );
is( $hsp->bits,                 '54.0' );
is( int $hsp->percent_identity, 88 );
float_is( $hsp->evalue, 2e-12 );
is( int( $hsp->frac_conserved * $hsp->length ), 83 );
is( $hsp->query->start,                         94 );
is( $hsp->query->end,                           180 );
is( $hsp->query->strand,                        1 );
is( $hsp->hit->strand,                          -1 );
is( $hsp->hit->start,                           1 );
is( $hsp->hit->end,                             94 );
is( $hsp->gaps,                                 7 );
is( $hsp->n,                                    1 );

# this is blastn bl2seq
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('bl2seq.blastn')
);
$result = $searchio->next_result;
isa_ok( $result, 'Bio::Search::Result::ResultI' );
is( $result->query_name,                          '' );
is( $result->query_length,                        180 );
is( $result->algorithm,                           'BLASTN' );
is( $result->algorithm_reference,                 undef );
is( $result->get_statistic('querylength'),        174 );
is( $result->get_statistic('effectivedblength'),  173 );
is( $result->get_statistic('effectivespace'),     30102 );
is( $result->get_statistic('effectivespaceused'), 30102 );
$hit = $result->next_hit;
is( $hit->name,   'human' );
is( $hit->length, 179 );
$hsp = $hit->next_hsp;
is( $hsp->score,                27 );
is( $hsp->bits,                 '54.0' );
is( int $hsp->percent_identity, 88 );
float_is( $hsp->evalue, 2e-12 );
is( int( $hsp->frac_conserved * $hsp->length ), 83 );
is( $hsp->query->start,                         94 );
is( $hsp->query->end,                           180 );
is( $hsp->query->strand,                        1 );
is( $hsp->hit->strand,                          1 );
is( $hsp->hit->start,                           86 );
is( $hsp->hit->end,                             179 );
is( $hsp->gaps,                                 7 );
is( $hsp->n,                                    1 );

# this is blastn bl2seq+
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('bl2seq+.blastn')
);
$result = $searchio->next_result;
isa_ok( $result, 'Bio::Search::Result::ResultI' );
is( $result->query_name, 'gi|2695846|emb|Y13255.1|' );
is( $result->query_description,
   'Acipenser baeri mRNA for immunoglobulin heavy chain, clone ScH 3.3'
);
is( $result->query_length,                        606 );
is( $result->algorithm,                          'BLASTN' );
is( $result->algorithm_version,                  '2.2.29+' );
is( $result->algorithm_reference,                 undef );
is( $result->get_statistic('effectivespaceused'), 352836 );
is( $result->get_statistic('kappa'),              0.621 );
is( $result->get_statistic('kappa_gapped'),      '0.460' );
is( $result->get_statistic('lambda'),             1.33 );
is( $result->get_statistic('lambda_gapped'),      1.28 );
is( $result->get_statistic('entropy'),            1.12 );
is( $result->get_statistic('entropy_gapped'),    '0.850' );
$hit = $result->next_hit;
is( $hit->name,   'gi|2695846|emb|Y13255.1|' );
is( $hit->description,
   'Acipenser baeri mRNA for immunoglobulin heavy chain, clone ScH 3.3'
);
is( $hit->length, 606 );
$hsp = $hit->next_hsp;
is( $hsp->score,            606 );
is( $hsp->bits,             1120 );
is( $hsp->percent_identity, 100 );
float_is( $hsp->evalue,    '0.0' );
is( $hsp->query->start,     1 );
is( $hsp->query->end,       606 );
is( $hsp->query->strand,    1 );
is( $hsp->hit->strand,      1 );
is( $hsp->hit->start,       1 );
is( $hsp->hit->end,         606 );
is( $hsp->gaps,             0 );
is( $hsp->n,                1 );

# this is blastp bl2seq
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('bl2seq.bug940.out')
);
$result = $searchio->next_result;
isa_ok( $result, 'Bio::Search::Result::ResultI' );
is( $result->query_name, 'zinc' );
is( $result->algorithm,  'BLASTP' );
is( $result->query_description,
'finger protein 135 (clone pHZ-17) [Homo sapiens]. neo_id RS.ctg14243-000000.6.0'
);
is( $result->query_length,                        469 );
is( $result->get_statistic('querylength'),        446 );
is( $result->get_statistic('effectivedblength'),  446 );
is( $result->get_statistic('effectivespace'),     198916 );
is( $result->get_statistic('effectivespaceused'), 198916 );
$hit = $result->next_hit;
is( $hit->name,    'gi|4507985|' );
is( $hit->ncbi_gi, 4507985 );
is( $hit->description,
'zinc finger protein 135 (clone pHZ-17) [Homo sapiens]. neo_id RS.ctg14243-000000.6.0'
);
is( $hit->length, 469 );
$hsp = $hit->next_hsp;
is( $hsp->score,                1626 );
is( $hsp->bits,                 637 );
is( int $hsp->percent_identity, 66 );
float_is( $hsp->evalue, 0.0 );
is( int( $hsp->frac_conserved * $hsp->length ), 330 );
is( $hsp->query->start,                         121 );
is( $hsp->query->end,                           469 );
is( $hsp->hit->start,                           1 );
is( $hsp->hit->end,                             469 );
is( $hsp->gaps,                                 120 );
is( $hsp->n,                                    1 );

ok( $hit->next_hsp );    # there is more than one HSP here,
                         # make sure it is parsed at least

# cannot distinguish between blastx and tblastn reports
# so we're only testing a blastx report for now

# this is blastx bl2seq
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('bl2seq.blastx.out')
);
$result = $searchio->next_result;
isa_ok( $result, 'Bio::Search::Result::ResultI' );
is( $result->query_name, 'AE000111.1' );
is( $result->query_description,
    'Escherichia coli K-12 MG1655 section 1 of 400 of the complete genome' );
is( $result->algorithm,                           'BLASTX' );
is( $result->algorithm_reference,                 undef );
is( $result->query_length,                        720 );
is( $result->get_statistic('querylength'),        undef );
is( $result->get_statistic('effectivedblength'),  787 );
is( $result->get_statistic('effectivespace'),     undef );
is( $result->get_statistic('effectivespaceused'), 162122 );
$hit = $result->next_hit;
is( $hit->name, 'AK1H_ECOLI' );
is( $hit->description,
'P00561 Bifunctional aspartokinase/homoserine dehydrogenase I (AKI-HDI) [Includes: Aspartokinase I ; Homoserine dehydrogenase I ]'
);
is( $hit->length, 820 );
$hsp = $hit->next_hsp;
is( $hsp->score,                634 );
is( $hsp->bits,                 248 );
is( int $hsp->percent_identity, 100 );
float_is( $hsp->evalue, 2e-70 );
is( int( $hsp->frac_conserved * $hsp->length ), 128 );
is( $hsp->query->start,                         1 );
is( $hsp->query->end,                           384 );
is( $hsp->hit->start,                           1 );
is( $hsp->hit->end,                             128 );
is( $hsp->gaps,                                 0 );
is( $hsp->query->frame,                         0 );
is( $hsp->hit->frame,                           0 );
is( $hsp->query->strand,                        -1 );
is( $hsp->hit->strand,                          0 );
is( $hsp->n,                                    1 );

# this is tblastx bl2seq (self against self)
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('bl2seq.tblastx.out')
);
$result = $searchio->next_result;
isa_ok( $result, 'Bio::Search::Result::ResultI' );
is( $result->query_name,          'Escherichia' );
is( $result->algorithm,           'TBLASTX' );
is( $result->algorithm_reference, undef );
is( $result->query_description,
    'coli K-12 MG1655 section 1 of 400 of the complete genome' );
is( $result->query_length,                        720 );
is( $result->get_statistic('querylength'),        undef );
is( $result->get_statistic('effectivedblength'),  221 );
is( $result->get_statistic('effectivespace'),     undef );
is( $result->get_statistic('effectivespaceused'), 48620 );
$hit = $result->next_hit;
is( $hit->name,    'gi|1786181|gb|AE000111.1|AE000111' );
is( $hit->ncbi_gi, 1786181 );
is( $hit->description,
    'Escherichia coli K-12 MG1655 section 1 of 400 of the complete genome' );
is( $hit->length, 720 );
$hsp = $hit->next_hsp;
is( $hsp->score,                1118 );
is( $hsp->bits,                 515 );
is( int $hsp->percent_identity, 95 );
float_is( $hsp->evalue, 1e-151 );
is( int( $hsp->frac_conserved * $hsp->length ), 229 );
is( $hsp->query->start,                         1 );
is( $hsp->query->end,                           720 );
is( $hsp->hit->start,                           1 );
is( $hsp->hit->end,                             720 );
is( $hsp->gaps,                                 0 );
is( $hsp->query->frame,                         0 );
is( $hsp->hit->frame,                           0 );
is( $hsp->query->strand,                        1 );
is( $hsp->hit->strand,                          1 );
is( $hsp->n,                                    1 );

# this is NCBI tblastn
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('tblastn.out')
);
$result = $searchio->next_result;
isa_ok( $result, 'Bio::Search::Result::ResultI' );
is( $result->algorithm, 'TBLASTN' );
like($result->algorithm_reference,qr/Gapped BLAST and PSI-BLAST/);
is( $result->get_statistic('querylength'),        102 );
is( $result->get_statistic('effectivedblength'),  4342 );
is( $result->get_statistic('effectivespace'),     442884 );
is( $result->get_statistic('effectivespaceused'), 442884 );
$hit = $result->next_hit;
is( $hit->name, 'gi|10040111|emb|AL390796.6|AL390796' );

# Test Blast parsing with B=0 (WU-BLAST)
$searchio = Bio::SearchIO->new(
    -file   => test_input_file('no_hsps.blastp'),
    -format => 'blast'
);
$result = $searchio->next_result;
like($result->algorithm_reference,qr/Gish, W. \(1996-2003\)/);
is( $result->query_name, 'mgri:MG00189.3' );
$hit = $result->next_hit;
is( $hit->name,        'mgri:MG00189.3' );
is( $hit->description, 'hypothetical protein 6892 8867 +' );
is( $hit->bits,        3098 );
float_is( $hit->significance, 0. );

$hit = $result->next_hit;
is( $hit->name,        'fgram:FG01141.1' );
is( $hit->description, 'hypothetical protein 47007 48803 -' );
is( $hit->bits,        2182 );
float_is( $hit->significance, 4.2e-226 );
is( $result->num_hits, 415 );

# Let's now test if _guess_format is doing its job correctly
my %pair = (
    'filename.blast' => 'blast',
    'filename.bls'   => 'blast',
    'f.blx'          => 'blast',
    'f.tblx'         => 'blast',
    'fast.bls'       => 'blast',
    'f.fasta'        => 'fasta',
    'f.fa'           => 'fasta',
    'f.fx'           => 'fasta',
    'f.fy'           => 'fasta',
    'f.ssearch'      => 'fasta',
    'f.SSEARCH.m9'   => 'fasta',
    'f.m9'           => 'fasta',
    'f.psearch'      => 'fasta',
    'f.osearch'      => 'fasta',
    'f.exon'         => 'exonerate',
    'f.exonerate'    => 'exonerate',
    'f.blastxml'     => 'blastxml',
    'f.xml'          => 'blastxml'
);
while ( my ( $file, $expformat ) = each %pair ) {
    is( Bio::SearchIO->_guess_format($file),
        $expformat, "$expformat for $file" );
}

# Test Wes Barris's reported bug when parsing blastcl3 output which
# has integer overflow

$searchio = Bio::SearchIO->new(
    -file   => test_input_file('hsinsulin.blastcl3.blastn'),
    -format => 'blast'
);
$result = $searchio->next_result;
is( $result->query_name,         'human' );
is( $result->database_letters(), '-24016349' );

# this is of course not the right length, but is the what blastcl3
# reports, the correct value is
is( $result->get_statistic('dbletters'), '192913178' );
is( $result->get_statistic('dbentries'), '1867771' );

# test for links and groups being parsed out of WU-BLAST properly
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('brassica_ATH.WUBLASTN')
);
ok( $result = $searchio->next_result );
ok( $hit    = $result->next_hit );
ok( $hsp    = $hit->next_hsp );
is( $hsp->links,         '(1)-3-2' );
is( $hsp->query->strand, 1 );
is( $hsp->hit->strand,   1 );
is( $hsp->hsp_group,     '1' );
is( $hsp->n,             1 );
## Web blast result parsing

$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('catalase-webblast.BLASTP')
);
ok( $result = $searchio->next_result );
is( $result->rid, '1118324516-16598-103707467515.BLASTQ1' );
ok( $hit = $result->next_hit );
is( $hit->name,      'gi|40747822|gb|EAA66978.1|', 'full hit name' );
is( $hit->accession, 'EAA66978',                   'hit accession' );
is( $hit->ncbi_gi,   40747822 );
ok( $hsp = $hit->next_hsp );
is( $hsp->query->start, 1,   'query start' );
is( $hsp->query->end,   528, 'query start' );
is( $hsp->n,            1 );

# tests for new BLAST 2.2.13 output
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('new_blastn.txt')
);

$result = $searchio->next_result;
is( $result->database_name,
'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS,GSS,environmental samples or phase 0, 1 or 2 HTGS sequences)'
);
is( $result->database_entries,  3742891 );
is( $result->database_letters,  16670205594 );
is( $result->algorithm,         'BLASTN' );
is( $result->algorithm_version, '2.2.13 [Nov-27-2005]' );
like($result->algorithm_reference, qr/Gapped BLAST and PSI-BLAST/);
is( $result->rid,                    '1141079027-8324-8848328247.BLASTQ4' );
is( $result->query_name,             'pyrR,' );
is( $result->query_length,           558 );
is( $result->get_statistic('kappa'), '0.711' );
is( $result->get_statistic('kappa_gapped'),        '0.711' );
is( $result->get_statistic('lambda'),              '1.37' );
is( $result->get_statistic('lambda_gapped'),       '1.37' );
is( $result->get_statistic('entropy'),             '1.31' );
is( $result->get_statistic('entropy_gapped'),      '1.31' );
is( $result->get_statistic('dbletters'),           '-509663586' );
is( $result->get_statistic('dbentries'),           3742891 );
is( $result->get_statistic('effective_hsplength'), undef );
is( $result->get_statistic('effectivespace'),      8935230198384 );
is(
    $result->get_statistic(
        'number_of_hsps_better_than_expect_value_cutoff_without_gapping'),
    0
);
is( $result->get_statistic('number_of_hsps_gapped'),              1771 );
is( $result->get_statistic('number_of_hsps_successfully_gapped'), 0 );
is( $result->get_statistic('length_adjustment'),                  22 );
is( $result->get_statistic('querylength'),                        536 );
is( $result->get_statistic('effectivedblength'),                  16670205594 );
is( $result->get_statistic('effectivespaceused'), 8891094027712 );
is( $result->get_parameter('matrix'),             'blastn matrix:1 -3' );
is( $result->get_parameter('gapopen'),            5 );
is( $result->get_parameter('gapext'),             2 );
is( $result->get_statistic('S2'),                 '60' );
is( $result->get_statistic('S2_bits'),            '119.4' );
float_is( $result->get_parameter('expect'), '1e-23' );
is( $result->get_statistic('num_extensions'), '117843' );

@valid = (
    [
        'gi|41400296|gb|AE016958.1|', 4829781, 'AE016958', 41400296, '6e-059',
        119, 236
    ],
    [
        'gi|54013472|dbj|AP006618.1|', 6021225, 'AP006618', 54013472, '4e-026',
        64, 127
    ],
    [
        'gi|57546753|dbj|BA000030.2|', 9025608, 'BA000030', 57546753, '1e-023',
        60, 119
    ]
);
$count = 0;

while ( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is( $hit->name,      shift @$d );
    is( $hit->length,    shift @$d );
    is( $hit->accession, shift @$d );
    is( $hit->ncbi_gi,   shift @$d );
    float_is( $hit->significance, shift @$d );
    is( $hit->raw_score, shift @$d );
    is( $hit->bits,      shift @$d );

    if ( $count == 0 ) {
        my $hsps_left = 1;
        while ( my $hsp = $hit->next_hsp ) {
            is( $hsp->query->start,    262 );
            is( $hsp->query->end,      552 );
            is( $hsp->hit->start,      1166897 );
            is( $hsp->hit->end,        1167187 );
            is( $hsp->length('total'), 291 );
            is( $hsp->hit_features,    'PyrR' );
            is( $hsp->start('hit'),    $hsp->hit->start );
            is( $hsp->end('query'),    $hsp->query->end );
            is( $hsp->strand('sbjct'), $hsp->subject->strand );  # alias for hit
            float_is( $hsp->evalue, 6e-59 );
            is( $hsp->score, 119 );
            is( $hsp->bits,  236 );
            is( sprintf( "%.2f", $hsp->percent_identity ),        85.22 );
            is( sprintf( "%.4f", $hsp->frac_identical('query') ), 0.8522 );
            is( sprintf( "%.4f", $hsp->frac_identical('hit') ),   0.8522 );
            is( $hsp->gaps, 0 );
            is( $hsp->n,    1 );
            $hsps_left--;
        }
        is( $hsps_left, 0 );
    }
    last if ( $count++ > @valid );
}
is( @valid, 0 );

# Bug 2189
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('blastp2215.blast')
);

$result = $searchio->next_result;
is( $result->database_entries,  4460989 );
is( $result->database_letters,  1533424333 );
is( $result->algorithm,         'BLASTP' );
is( $result->algorithm_version, '2.2.15 [Oct-15-2006]' );
is( $result->rid,               '1169055516-21385-22799250964.BLASTQ4' );
is( $result->query_name,        'gi|15608519|ref|NP_215895.1|' );
is( $result->query_gi,          15608519 );
is( $result->query_length,      193 );
@hits = $result->hits;
is( scalar(@hits),          10 );
is( $hits[1]->accession,    '1W30' );
is( $hits[4]->significance, '2e-72' );
is( $hits[7]->bits,         '254' );
$result = $searchio->next_result;
is( $result->database_entries,  4460989 );
is( $result->database_letters,  1533424333 );
is( $result->algorithm,         'BLASTP' );
is( $result->algorithm_version, '2.2.15 [Oct-15-2006]' );
is( $result->query_name,        'gi|15595598|ref|NP_249092.1|' );
is( $result->query_length,      423 );
@hits = $result->hits;
is( scalar(@hits),          10 );
is( $hits[1]->accession,    'ZP_00972546' );
is( $hits[2]->ncbi_gi,      116054132 );
is( $hits[4]->significance, '0.0' );
is( $hits[7]->bits,         624 );

# Bug 2246
$searchio = Bio::SearchIO->new(
    -format  => 'blast',
    -verbose => -1,
    -file    => test_input_file('bug2246.blast')
);
$result = $searchio->next_result;
is(
    $result->get_statistic(
        'number_of_hsps_better_than_expect_value_cutoff_without_gapping'),
    0
);
is( $result->get_statistic('number_of_hsps_gapped'),              7049 );
is( $result->get_statistic('number_of_hsps_successfully_gapped'), 55 );
is( $result->get_statistic('length_adjustment'),                  125 );
is( $result->get_statistic('querylength'),                        68 );
is( $result->get_statistic('effectivedblength'),                  1045382588 );
is( $result->get_statistic('effectivespace'),                     71086015984 );
is( $result->get_statistic('effectivespaceused'),                 71086015984 );
$hit = $result->next_hit;
is $hit->name,        'UniRef50_Q9X0H5';
is $hit->length,      0;
is $hit->accession,   'UniRef50_Q9X0H5';
is $hit->description, 'Cluster: Histidyl-tRNA synthetase; n=4; Thermoto...';
is $hit->bits,        23;
float_is( $hit->significance, 650 );

# Bug 1986
$searchio = Bio::SearchIO->new(
    -format  => 'blast',
    -verbose => -1,
    -file    => test_input_file('bug1986.blastp')
);
$result = $searchio->next_result;
is( $result->get_statistic('querylength'),        335 );
is( $result->get_statistic('effectivedblength'),  18683311 );
is( $result->get_statistic('effectivespace'),     6258909185 );
is( $result->get_statistic('effectivespaceused'), 6258909185 );
$hit = $result->next_hit;
is $hit->name,      'ENSP00000350182';
is $hit->length,    425;
is $hit->accession, 'ENSP00000350182';
is $hit->description,
'pep:novel clone::BX322644.8:4905:15090:-1 gene:ENSG00000137397 transcript:ENST00000357569';
is $hit->raw_score, 301;
is $hit->bits,      120;
float_is( $hit->significance, 3e-27 );
$hit = $result->next_hit;
is $hit->name,      'ENSP00000327738';
is $hit->length,    468;
is $hit->accession, 'ENSP00000327738';
is $hit->description,
'pep:known-ccds chromosome:NCBI36:4:189297592:189305643:1 gene:ENSG00000184108 transcript:ENST00000332517 CCDS3851.1';
is $hit->raw_score, 289;
is $hit->bits,      115;
float_is( $hit->significance, 8e-26 );

# Bug 1986, pt. 2

# handle at least the first iteration with BLAST searches using databases
# containing non-unique IDs

my $file = test_input_file('bug1986.blast2');
my %unique_accs;
open my $IN, '<', $file or die "Could not read file '$file': $!\n";

while (<$IN>) {
    last if (/^Sequences/);
}
$count = 1;
while (<$IN>) {
    chomp;
    next if m{^\s*$};
    next unless ($_);
    last if m{^>};
    my ($accession) = split(/\s+/);

    #print "Real Hit $count = $accession\n";
    $unique_accs{$accession}++;

    #last if ($count == 10);
    ++$count;
}
close $IN;

is( $count,                      495 );
is( scalar( keys %unique_accs ), 490 );

my %search_accs;

$searchio = Bio::SearchIO->new(
    -format  => 'blast',
    -verbose => -1,
    -file    => $file
);
$result = $searchio->next_result;
$count  = 1;
while ( my $hit = $result->next_hit ) {
    $search_accs{ $hit->accession }++;
    $count++;
}

is( $count,                      495 );
is( scalar( keys %search_accs ), 490 );

is_deeply( \%unique_accs, \%search_accs );

# bug 2391 - long query names

$file = test_input_file('bug2391.megablast');

$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => $file
);
$result = $searchio->next_result;

# data is getting munged up with long names
is( $result->query_name,
    'c6_COX;c6_QBL;6|31508172;31503325;31478402|rs36223351|1|dbSNP|C/G' );
is( $result->query_description, '' );
is( $result->algorithm,         'MEGABLAST' );
is(
    $result->get_statistic(
        'number_of_hsps_better_than_expect_value_cutoff_without_gapping'),
    undef
);
is( $result->get_statistic('number_of_hsps_gapped'),              0 );
is( $result->get_statistic('number_of_hsps_successfully_gapped'), 0 );
is( $result->get_statistic('length_adjustment'),                  16 );
is( $result->get_statistic('querylength'),                        85 );
is( $result->get_statistic('effectivedblength'),                  59358266 );
is( $result->get_statistic('effectivespace'),                     5045452610 );
is( $result->get_statistic('effectivespaceused'),                 5045452610 );

# bug 2399 - catching Expect(n) values

$file = test_input_file('bug2399.tblastn');

$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => $file
);
my $total_n = 0;
while ( my $query = $searchio->next_result ) {
    while ( my $subject = $query->next_hit ) {
        $total_n += grep { $_->n } $subject->hsps;
    }
}
is( $total_n, 80 );    # n = at least 1, so this was changed to reflect that

sub cmp_evalue ($$) {
    my ( $tval, $aval ) = @_;
    is( sprintf( "%g", $tval ), sprintf( "%g", $aval ) );
}

# bug 3064 - All-gap Query/Subject lines for BLAST+ do not have numbering

$file = test_input_file('blast_plus.blastp');

$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => $file
);

my $total_hsps = 0;
while ( my $query = $searchio->next_result ) {
    is( $query->get_statistic('querylength'),        undef );
    is( $query->get_statistic('effectivedblength'),  undef );
    is( $query->get_statistic('effectivespace'),     undef );
    is( $query->get_statistic('effectivespaceused'), 55770 );
    while ( my $subject = $query->next_hit ) {
        while ( my $hsp = $subject->next_hsp ) {
            $total_hsps++;
            if ( $total_hsps == 1 ) {
                is( $hsp->start('query'),         5 );
                is( $hsp->start('hit'),           3 );
                is( $hsp->end('query'),           220 );
                is( $hsp->end('hit'),             308 );
                is( length( $hsp->query_string ), length( $hsp->hit_string ) );
            }
        }
    }
}

is( $total_hsps, 2 );

# BLAST 2.2.20+ output file ZABJ4EA7014.CH878695.1.blast.txt
# Tests SearchIO blast parsing of 'Features in/flanking this part of a subject sequence'
$searchio = Bio::SearchIO->new(
    -format => 'blast',
    -file   => test_input_file('ZABJ4EA7014.CH878695.1.blast.txt')
);

$result = $searchio->next_result;

# Parse BLAST header details
is( $result->algorithm,         'BLASTN' );
is( $result->algorithm_version, '2.2.20+' );
like($result->algorithm_reference, qr/A greedy algorithm for aligning DNA\s+sequences/);
is( $result->database_name,
    'human build 35 genome database (reference assembly only)' );
is( $result->database_entries, 378 );
is( $result->database_letters, 2866055344 );
is( $result->query_name,       'gi|95131563|gb|CH878695.1|' );
is( $result->query_description,
    'Homo sapiens 211000035829648 genomic scaffold' );
is( $result->query_length, 29324 );

# Parse BLAST footer details
is( $result->get_statistic('posted_date'),    'Jul 26, 2007  3:20 PM' );
is( $result->get_statistic('dbletters'),      -1428911948 );
is( $result->get_statistic('lambda'),         '1.33' );
is( $result->get_statistic('kappa'),          '0.621' );
is( $result->get_statistic('entropy'),        '1.12' );
is( $result->get_statistic('lambda_gapped'),  '1.28' );
is( $result->get_statistic('kappa_gapped'),   '0.460' );
is( $result->get_statistic('entropy_gapped'), '0.850' );
is( $result->get_parameter('matrix'),         'blastn matrix:1 -2' );
is( $result->get_parameter('gapopen'),        0 );
is( $result->get_parameter('gapext'),         0 );
is( $result->get_statistic('num_extensions'), 216 );
is( $result->get_statistic('num_successful_extensions'), 216 );
is( $result->get_parameter('expect'),                    '0.01' );
is( $result->get_statistic('seqs_better_than_cutoff'),   10 );
is(
    $result->get_statistic(
        'number_of_hsps_better_than_expect_value_cutoff_without_gapping'),
    0
);
is( $result->get_statistic('number_of_hsps_gapped'),              216 );
is( $result->get_statistic('number_of_hsps_successfully_gapped'), 212 );
is( $result->get_statistic('length_adjustment'),                  34 );
is( $result->get_statistic('querylength'),                        29290 );
is( $result->get_statistic('effectivedblength'),                  2866042492 );
is( $result->get_statistic('effectivespace'),     83946384590680 );
is( $result->get_statistic('effectivespaceused'), 83946384590680 );
is( $result->get_statistic('A'),                  0 );
is( $result->get_statistic('X1'),                 23 );
is( $result->get_statistic('X1_bits'),            '44.2' );
is( $result->get_statistic('X2'),                 32 );
is( $result->get_statistic('X2_bits'),            '59.1' );
is( $result->get_statistic('X3'),                 54 );
is( $result->get_statistic('X3_bits'),            '99.7' );
is( $result->get_statistic('S1'),                 23 );
is( $result->get_statistic('S1_bits'),            '43.6' );
is( $result->get_statistic('S2'),                 29 );
is( $result->get_statistic('S2_bits'),            '54.7' );

# Skip the 1st hit. It doesn't have any 'Features in/flanking this part of subject sequence:'
$hit = $result->next_hit;

# The 2nd hit has hsps with 'Features flanking this part of subject sequence:'
$hit = $result->next_hit;
is( $hit->name,        'gi|51459264|ref|NT_077382.3|Hs1_77431' );
is( $hit->description, 'Homo sapiens chromosome 1 genomic contig' );
is( $hit->length,      237250 );

# In the 2nd hit, look at the 1st hsp
$hsp = $hit->next_hsp;
is( $hsp->hit_features,
"16338 bp at 5' side: PRAME family member 8   11926 bp at 3' side: PRAME family member 9"
);
is( $hsp->bits,          7286 );
is( $hsp->score,         3945 );
is( $hsp->expect,        '0.0' );
is( $hsp->hsp_length,    6145 );
is( $hsp->num_identical, 5437 );
is( int sprintf( "%.2f", $hsp->percent_identity ), 88 );
is( $hsp->gaps,           152 );
is( $hsp->start('query'), 23225 );
is( $hsp->start('sbjct'), 86128 );
is( $hsp->end('query'),   29324 );
is( $hsp->end('sbjct'),   92165 );

# In the 2nd hit, look at the 2nd hsp
$hsp = $hit->next_hsp;
is( $hsp->hit_features,
"25773 bp at 5' side: PRAME family member 3   3198 bp at 3' side: PRAME family member 5"
);
is( $hsp->bits,          4732 );
is( $hsp->score,         2562 );
is( $hsp->expect,        '0.0' );
is( $hsp->hsp_length,    4367 );
is( $hsp->num_identical, 3795 );
is( int sprintf( "%.2f", $hsp->percent_identity ), 86 );
is( $hsp->gaps,           178 );
is( $hsp->start('query'), 23894 );
is( $hsp->start('sbjct'), 37526 );
is( $hsp->end('query'),   28193 );
is( $hsp->end('sbjct'),   41781 );

# In the 2nd hit, look at the 3rd hsp
$hsp = $hit->next_hsp;
is( $hsp->hit_features,
"16338 bp at 5' side: PRAME family member 8   14600 bp at 3' side: PRAME family member 9"
);
is( $hsp->bits,          3825 );
is( $hsp->score,         2071 );
is( $hsp->expect,        '0.0' );
is( $hsp->hsp_length,    3406 );
is( $hsp->num_identical, 2976 );
is( int sprintf( "%.2f", $hsp->percent_identity ), 87 );
is( $hsp->gaps,           89 );
is( $hsp->start('query'), 14528 );
is( $hsp->start('sbjct'), 86128 );
is( $hsp->end('query'),   17886 );
is( $hsp->end('sbjct'),   89491 );

# In the 2nd hit, look at the 4th hsp
$hsp = $hit->next_hsp;
is( $hsp->hit_features,
"29101 bp at 5' side: PRAME family member 8   2120 bp at 3' side: PRAME family member 9"
);
is( $hsp->bits,          3241 );
is( $hsp->score,         1755 );
is( $hsp->expect,        '0.0' );
is( $hsp->hsp_length,    3158 );
is( $hsp->num_identical, 2711 );
is( int sprintf( "%.2f", $hsp->percent_identity ), 85 );
is( $hsp->gaps,           123 );
is( $hsp->start('query'), 23894 );
is( $hsp->start('sbjct'), 98891 );
is( $hsp->end('query'),   27005 );
is( $hsp->end('sbjct'),   101971 );

# In the 2nd hit, look at the 5th hsp
$hsp = $hit->next_hsp;
is( $hsp->hit_features,  "PRAME family member 13" );
is( $hsp->bits,          3142 );
is( $hsp->score,         1701 );
is( $hsp->expect,        '0.0' );
is( $hsp->hsp_length,    2507 );
is( $hsp->num_identical, 2249 );
is( int sprintf( "%.2f", $hsp->percent_identity ), 89 );
is( $hsp->gaps,           63 );
is( $hsp->start('query'), 3255 );
is( $hsp->start('sbjct'), 128516 );
is( $hsp->end('query'),   5720 );
is( $hsp->end('sbjct'),   131000 );

# testing for Bug #3298
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('multiresult_blastn+.bls')
);

is ($searchio->next_result->algorithm_version, '2.2.25+', "testing Bug 3298");
is ($searchio->next_result->algorithm_version, '2.2.25+', "testing Bug 3298");
is ($searchio->next_result->algorithm_version, '2.2.25+', "testing Bug 3298");

# testing for Bug #3251
$searchio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('rpsblast_no_hits.bls')
);

is ($searchio->next_result->database_name, 'CDD.v.2.13', "testing Bug 3251");
is ($searchio->next_result->database_name, 'CDD.v.2.13', "testing Bug 3251");
is ($searchio->next_result->database_name, 'CDD.v.2.13', "testing Bug 3251");
