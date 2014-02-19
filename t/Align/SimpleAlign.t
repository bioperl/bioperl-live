# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 206 );

    use_ok('Bio::SimpleAlign');
    use_ok('Bio::AlignIO');
    use_ok('Bio::SeqFeature::Generic');
    use_ok('Bio::Location::Simple');
    use_ok('Bio::Location::Split');
}

my $DEBUG = test_debug();

my ( $str, $aln, @seqs, $seq );

$str = Bio::AlignIO->new( -file => test_input_file('testaln.pfam') );
isa_ok( $str, 'Bio::AlignIO' );
$aln = $str->next_aln();
is $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', "pfam input test";

my $aln1 = $aln->remove_columns( ['mismatch'] );
is(
    $aln1->match_line,
    '::*::::*:**:*:*:***:**.***::*.*::**::**:***..**:'
      . '*:*.::::*:.:*.*.**:***.**:*.:.**::**.*:***********:::*:.:*:**.*::*:'
      . '.*.:*:**:****************::',
    'match_line'
);

my $aln2 = $aln->select( 1, 3 );
isa_ok( $aln2, 'Bio::Align::AlignI' );
is( $aln2->num_sequences, 3, 'num_sequences' );

# test select non contiguous-sorted by default
$aln2 = $aln->select_noncont( 6, 7, 8, 9, 10, 1, 2, 3, 4, 5 );
is( $aln2->num_sequences, 10, 'num_sequences' );
is(
    $aln2->get_seq_by_pos(2)->id,
    $aln->get_seq_by_pos(2)->id,
    'select_noncont'
);
is(
    $aln2->get_seq_by_pos(8)->id,
    $aln->get_seq_by_pos(8)->id,
    'select_noncont'
);

# test select non contiguous-nosort option
$aln2 = $aln->select_noncont( 'nosort', 6, 7, 8, 9, 10, 1, 2, 3, 4, 5 );
is( $aln2->num_sequences, 10, 'num_sequences' );
is(
    $aln2->get_seq_by_pos(2)->id,
    $aln->get_seq_by_pos(7)->id,
    'select_noncont'
);
is(
    $aln2->get_seq_by_pos(8)->id,
    $aln->get_seq_by_pos(3)->id,
    'select_noncont'
);

# test select non contiguous by name
my $aln3 = $aln->select_noncont_by_name('1433_LYCES','BMH1_YEAST','143T_HUMAN');
is( $aln3->num_sequences, 3, 'select_noncont_by_name' );
my @seqs3 = $aln3->each_seq();
is $seqs3[0]->id, '1433_LYCES', 'select_noncont_by_name';
is $seqs3[1]->id, 'BMH1_YEAST', 'select_noncont_by_name';
is $seqs3[2]->id, '143T_HUMAN', 'select_noncont_by_name';


@seqs = $aln->each_seq();
is scalar @seqs, 16, 'each_seq';
is $seqs[0]->get_nse, '1433_LYCES/9-246', 'get_nse';
is $seqs[0]->id,      '1433_LYCES',       'id';
is $seqs[0]->num_gaps, 3,                  'num_gaps';
@seqs = $aln->each_alphabetically();
is scalar @seqs, 16, 'each_alphabetically';

is $aln->column_from_residue_number( '1433_LYCES', 10 ), 2,
  'column_from_residue_number';
is $aln->displayname( '1433_LYCES/9-246', 'my_seq' ), 'my_seq',
  'display_name get/set';
is $aln->displayname('1433_LYCES/9-246'), 'my_seq', 'display_name get';
is substr( $aln->consensus_string(50), 0, 60 ),
  "RE??VY?AKLAEQAERYEEMV??MK?VAE??????ELSVEERNLLSVAYKNVIGARRASW",
  'consensus_string';
is substr( $aln->consensus_string(100), 0, 60 ),
  "?????????L????E????M???M????????????L??E?RNL?SV?YKN??G??R??W",
  'consensus_string';
is substr( $aln->consensus_string(0), 0, 60 ),
  "REDLVYLAKLAEQAERYEEMVEFMKKVAELGAPAEELSVEERNLLSVAYKNVIGARRASW",
  'consensus_string';

ok( @seqs = $aln->each_seq_with_id('143T_HUMAN') );
is scalar @seqs, 1, 'each_seq_with_id';

is $aln->is_flush, 1, 'is_flush';
ok( $aln->id('x') && $aln->id eq 'x', 'id get/set' );

is $aln->length,        242,  'length';
is $aln->num_residues,  3769, 'num_residues';
is $aln->num_sequences, 16,   'num_sequences';
is( sprintf( "%.2f", $aln->overall_percentage_identity() ),
    33.06, 'overall_percentage_identity' );
is( sprintf( "%.2f", $aln->overall_percentage_identity('align') ),
    33.06, 'overall_percentage_identity (align)' );
is( sprintf( "%.2f", $aln->overall_percentage_identity('short') ),
    35.24, 'overall_percentage_identity (short)' );
is( sprintf( "%.2f", $aln->overall_percentage_identity('long') ),
    33.47, 'overall_percentage_identity (long)' );
is( sprintf( "%.2f", $aln->average_percentage_identity() ),
    66.91, 'average_percentage_identity' );

ok $aln->set_displayname_count;
is $aln->displayname('1433_LYCES/9-246'), '1433_LYCES_1',
  'set_displayname_count';
ok $aln->set_displayname_flat;
is $aln->displayname('1433_LYCES/9-246'), '1433_LYCES', 'set_displayname_flat';
ok $aln->set_displayname_normal;
is $aln->displayname('1433_LYCES/9-246'), '1433_LYCES/9-246',
  'set_displayname_normal';
ok $aln->uppercase;
ok $aln->map_chars( '\.', '-' );
@seqs = $aln->each_seq_with_id('143T_HUMAN');
is substr( $seqs[0]->seq, 0, 60 ),
  'KTELIQKAKLAEQAERYDDMATCMKAVTEQGA---ELSNEERNLLSVAYKNVVGGRRSAW',
  'uppercase, map_chars';

is(
    $aln->match_line,
    '       ::*::::*  : *   *:           *: *:***:**.***::*.'
      . ' *::**::**:***      .  .      **  :* :*   .  :: ::   *:  .     :* .*. **:'
      . '***.** :*.            :  .*  *   :   : **.*:***********:::* : .: *  :** .'
      . '*::*: .*. : *: **:****************::     ',
    'match_line'
);
ok $aln->remove_seq( $seqs[0] ), 'remove_seqs';
is $aln->num_sequences, 15, 'remove_seqs';
ok $aln->add_seq( $seqs[0] ), 'add_seq';
is $aln->num_sequences, 16, 'add_seq';
ok $seq = $aln->get_seq_by_pos(1), 'get_seq_by_pos';
is( $seq->id, '1433_LYCES', 'get_seq_by_pos' );
ok( ( $aln->missing_char(), 'P' ) and ( $aln->missing_char('X'), 'X' ) );
ok( ( $aln->match_char(),   '.' ) and ( $aln->match_char('-'),   '-' ) );
ok( ( $aln->gap_char(),     '-' ) and ( $aln->gap_char('.'),     '.' ) );

is $aln->purge(0.7), 12, 'purge';
is $aln->num_sequences, 4, 'purge';

SKIP: {
    test_skip( -tests => 24, -requires_module => 'IO::String' );

    my $string;
    my $out = IO::String->new($string);

    my $s1 = Bio::LocatableSeq->new(
        -id       => 'AAA',
        -seq      => 'aawtat-tn-',
        -start    => 1,
        -end      => 8,
        -alphabet => 'dna'
    );
    my $s2 = Bio::LocatableSeq->new(
        -id       => 'BBB',
        -seq      => '-aaaat-tt-',
        -start    => 1,
        -end      => 7,
        -alphabet => 'dna'
    );
    $a = Bio::SimpleAlign->new();
    $a->add_seq($s1);
    $a->add_seq($s2);

    is( $a->consensus_iupac, "aAWWAT-TN-", 'IO::String consensus_iupac' );
    $s1->seq('aaaaattttt');
    $s1->alphabet('dna');
    $s1->end(10);
    $s2->seq('-aaaatttt-');
    $s2->end(8);
    $a = Bio::SimpleAlign->new();
    $a->add_seq($s1);
    $a->add_seq($s2);

    my $strout = Bio::AlignIO->new( -fh => $out, -format => 'pfam' );
    $strout->write_aln($a);
    is(
        $string,
        "AAA/1-10    aaaaattttt\n" . "BBB/1-8     -aaaatttt-\n",
        'IO::String write_aln normal'
    );

    $out->setpos(0);
    $string = '';
    my $b = $a->slice( 2, 9 );
    $strout->write_aln($b);
    is $string,
      "AAA/2-9    aaaatttt\n" . "BBB/1-8    aaaatttt\n",
      'IO::String write_aln slice';

    $out->setpos(0);
    $string = '';
    $b = $a->slice( 9, 10 );
    $strout->write_aln($b);
    is $string,
      "AAA/9-10    tt\n" . "BBB/8-8     t-\n",
      'IO::String write_aln slice';

    $a->verbose(-1);
    $out->setpos(0);
    $string = '';
    $b = $a->slice( 1, 2 );
    $strout->write_aln($b);
    is $string,
      "AAA/1-2    aa\n" . "BBB/1-1    -a\n",
      'IO::String write_aln slice';

    # not sure what coordinates this should return...
    $a->verbose(-1);
    $out->setpos(0);
    $string = '';
    $b = $a->slice( 1, 1, 1 );
    $strout->write_aln($b);
    is $string,
      "AAA/1-1    a\n" . "BBB/1-0    -\n",
      'IO::String write_aln slice';

    $a->verbose(-1);
    $out->setpos(0);
    $string = '';
    $b = $a->slice( 2, 2 );
    $strout->write_aln($b);
    is $string,
      "AAA/2-2    a\n" . "BBB/1-1    a\n",
      'IO::String write_aln slice';

    eval { $b = $a->slice( 11, 13 ); };

    like( $@, qr/EX/ );

    # remove_columns by position
    $out->setpos(0);
    $string = '';
    $str    = Bio::AlignIO->new( -file => test_input_file('mini-align.aln') );
    $aln1   = $str->next_aln;
    $aln2   = $aln1->remove_columns( [ 0, 0 ] );
    $strout->write_aln($aln2);
    is $string,
        "P84139/2-33              NEGEHQIKLDELFEKLLRARLIFKNKDVLRRC\n"
      . "P814153/2-33             NEGMHQIKLDVLFEKLLRARLIFKNKDVLRRC\n"
      . "BAB68554/1-14            ------------------AMLIFKDKQLLQQC\n"
      . "gb|443893|124775/1-32    MRFRFQIKVPPAVEGARPALLIFKSRPELGGC\n",
      'remove_columns by position';

    # and when arguments are entered in "wrong order"?
    $out->setpos(0);
    $string = '';
    my $aln3 = $aln1->remove_columns( [ 1, 1 ], [ 30, 30 ], [ 5, 6 ] );
    $strout->write_aln($aln3);
    is $string,
        "P84139/1-33              MEGEIKLDELFEKLLRARLIFKNKDVLRC\n"
      . "P814153/1-33             MEGMIKLDVLFEKLLRARLIFKNKDVLRC\n"
      . "BAB68554/1-14            ----------------AMLIFKDKQLLQC\n"
      . "gb|443893|124775/2-32    -RFRIKVPPAVEGARPALLIFKSRPELGC\n",
      'remove_columns by position (wrong order)';

    my %cigars = $aln1->cigar_line;
    is $cigars{'gb|443893|124775/1-32'}, '19,19:21,24:29,29:32,32',
      'cigar_line';
    is $cigars{'P814153/1-33'},  '20,20:22,25:30,30:33,33', 'cigar_line';
    is $cigars{'BAB68554/1-14'}, '1,1:3,6:11,11:14,14',     'cigar_line';
    is $cigars{'P84139/1-33'},   '20,20:22,25:30,30:33,33', 'cigar_line';

    # sort_alphabetically
    my $s3 = Bio::LocatableSeq->new(
        -id       => 'ABB',
        -seq      => '-attat-tt-',
        -start    => 1,
        -end      => 7,
        -alphabet => 'dna'
    );
    $a->add_seq($s3);

    is $a->get_seq_by_pos(2)->id, "BBB", 'sort_alphabetically - before';
    ok $a->sort_alphabetically;
    is $a->get_seq_by_pos(2)->id, "ABB", 'sort_alphabetically - after';

    $b = $a->remove_gaps();
    is $b->consensus_string, "aaaattt", 'remove_gaps';

    $s1->seq('aaaaattt--');

    $b = $a->remove_gaps( undef, 'all_gaps_only' );
    is $b->consensus_string, "aaaaatttt", 'remove_gaps all_gaps_only';

    # test set_new_reference:
    $str = Bio::AlignIO->new( -file => test_input_file('testaln.clustalw') );
    $aln = $str->next_aln();
    my $new_aln = $aln->set_new_reference(3);
    $a       = $new_aln->get_seq_by_pos(1)->display_id;
    $new_aln = $aln->set_new_reference('P851414');
    $b       = $new_aln->get_seq_by_pos(1)->display_id;
    is $a, 'P851414', 'set_new_reference';
    is $b, 'P851414', 'set_new_reference';

    # test uniq_seq:
    $str = Bio::AlignIO->new(
        -verbose => $DEBUG,
        -file    => test_input_file('testaln2.fasta')
    );
    $aln     = $str->next_aln();
    $new_aln = $aln->uniq_seq();
    $a       = $new_aln->num_sequences;
    is $a, 11, 'uniq_seq';

    # check if slice works well with a LocateableSeq in its negative strand
    my $seq1 = Bio::LocatableSeq->new(
        -SEQ    => "ATGCTG-ATG",
        -START  => 1,
        -END    => 9,
        -ID     => "test1",
        -STRAND => -1
    );

    my $seq2 = Bio::LocatableSeq->new(
        -SEQ    => "A-GCTGCATG",
        -START  => 1,
        -END    => 9,
        -ID     => "test2",
        -STRAND => 1
    );

    $string = '';
    my $aln_negative = Bio::SimpleAlign->new();
    $aln_negative->add_seq($seq1);
    $aln_negative->add_seq($seq2);
    my $start_column =
      $aln_negative->column_from_residue_number(
        $aln_negative->get_seq_by_pos(1)->display_id, 2 );
    my $end_column =
      $aln_negative->column_from_residue_number(
        $aln_negative->get_seq_by_pos(1)->display_id, 5 );
    $aln_negative = $aln_negative->slice( $end_column, $start_column );
    my $seq_negative = $aln_negative->get_seq_by_pos(1);
    is( $seq_negative->start, 2, "bug 2099" );
    is( $seq_negative->end,   5, "bug 2099" );

    # bug 2793
    my $s11 = Bio::LocatableSeq->new( -id => 'testseq1', -seq => 'AAA' );
    my $s21 = Bio::LocatableSeq->new( -id => 'testseq2', -seq => 'CCC' );
    $a = Bio::SimpleAlign->new();
    ok( $a->add_seq( $s11, 1 ), "bug 2793" );
    is( $a->get_seq_by_pos(1)->seq, 'AAA', "bug 2793" );
    ok( $a->add_seq( $s21, 2 ), "bug 2793" );
    is( $a->get_seq_by_pos(2)->seq, 'CCC', "bug 2793" );
    throws_ok { $a->add_seq( $s21, 0 ) } qr/must be >= 1/, 'Bad sequence, bad!';
}

# test for Bio::SimpleAlign annotation method and
# Bio::FeatureHolder stuff

$aln = Bio::SimpleAlign->new;
isa_ok($aln,"Bio::AnnotatableI");

for my $seqset ( [qw(one AGAGGAT)], [qw(two AGACGAT)], [qw(three AGAGGTT)] ) {
    $aln->add_seq(
        Bio::LocatableSeq->new(
            -id  => $seqset->[0],
            -seq => $seqset->[1]
        )
    );
}

is $aln->num_sequences, 3, 'added 3 seqs';

$aln->add_SeqFeature(
    Bio::SeqFeature::Generic->new(
        -start       => 1,
        -end         => 1,
        -primary_tag => 'charLabel',
    )
);
$aln->add_SeqFeature(
    Bio::SeqFeature::Generic->new(
        -start       => 3,
        -end         => 3,
        -primary_tag => 'charLabel',

    )
);
is( $aln->feature_count, 2, 'first 2 features added' );

my $splitloc = Bio::Location::Split->new;
$splitloc->add_sub_Location(
    Bio::Location::Simple->new(
        -start => 2,
        -end   => 3
    )
);

$splitloc->add_sub_Location(
    Bio::Location::Simple->new(
        -start => 5,
        -end   => 6
    )
);

$aln->add_SeqFeature(
    Bio::SeqFeature::Generic->new(
        -location    => $splitloc,
        -primary_tag => 'charLabel',
    )
);

is( $aln->feature_count, 3, '3rd feature added' );

#do slices and masks as defined by the feature
my $i          = 0;
my @slice_lens = qw(1 1 2 2);
for my $feature ( $aln->get_SeqFeatures ) {
    for my $loc ( $feature->location->each_Location ) {
        my $masked = $aln->mask_columns( $loc->start, $loc->end, '?');
        $masked->verbose(2);
        lives_ok {my $fslice = $masked->slice( $loc->start, $loc->end )};

        $masked->verbose(-1);
        my $fslice = $masked->slice( $loc->start, $loc->end );
        is( $fslice->length, $slice_lens[ $i++ ], "slice $i len" );
        for my $s ( $fslice->each_seq ) {
            like( $s->seq, qr/^\?+$/, 'correct masked seq' );
        }
    }
}

# test set_displayname_safe & restore_displayname:
$str = Bio::AlignIO->new( -file => test_input_file('pep-266.aln') );
$aln = $str->next_aln();
is $aln->get_seq_by_pos(3)->display_id, 'Smik_Contig1103.1',
  'initial display id ok';
my ( $new_aln, $ref ) = $aln->set_displayname_safe();
is $new_aln->get_seq_by_pos(3)->display_id, 'S000000003', 'safe display id ok';
my $restored_aln = $new_aln->restore_displayname($ref);
is $restored_aln->get_seq_by_pos(3)->display_id, 'Smik_Contig1103.1',
  'restored display id ok';

# test sort_by_list:
$str = Bio::AlignIO->new( -file => test_input_file('testaln.clustalw') );
my $list_file = test_input_file('testaln.list');
$aln     = $str->next_aln();
$new_aln = $aln->sort_by_list($list_file);
$a       = $new_aln->get_seq_by_pos(1)->display_id;
is $a, 'BAB68554', 'sort by list ok';

# test for Binary/Morphological/Mixed data

# sort_by_start

# test sort_by_list:

my $s1 = Bio::LocatableSeq->new(
    -id       => 'AAA',
    -seq      => 'aawtat-tn-',
    -start    => 12,
    -end      => 19,
    -alphabet => 'dna'
);
my $s2 = Bio::LocatableSeq->new(
    -id       => 'BBB',
    -seq      => '-aaaat-tt-',
    -start    => 1,
    -end      => 7,
    -alphabet => 'dna'
);
my $s3 = Bio::LocatableSeq->new(
    -id       => 'BBB',
    -seq      => '-aaaat-tt-',
    -start    => 31,
    -end      => 37,
    -alphabet => 'dna'
);
$a = Bio::SimpleAlign->new();
$a->add_seq($s1);
$a->add_seq($s2);
$a->add_seq($s3);

@seqs = $a->each_seq;
is( $seqs[0]->start, 12 );
is( $seqs[1]->start, 1 );
is( $seqs[2]->start, 31 );

$a->sort_by_start;
@seqs = $a->each_seq;

is( $seqs[0]->start, 1 );
is( $seqs[1]->start, 12 );
is( $seqs[2]->start, 31 );

my %testdata = (
    'allele1' => 'GGATCCATT[C/C]CTACT',
    'allele2' => 'GGAT[C/-][C/-]ATT[C/C]CT[A/C]CT',
    'allele3' => 'G[G/C]ATCCATT[C/G]CTACT',
    'allele4' => 'GGATCCATT[C/G]CTACT',
    'allele5' => 'GGATCCATT[C/G]CTAC[T/A]',
    'allele6' => 'GGATCCATT[C/G]CTA[C/G][T/A]',
    'testseq' => 'GGATCCATT[C/G]CTACT'
);

my $alnin = Bio::AlignIO->new(
    -format => 'fasta',
    -file   => test_input_file('alleles.fas')
);

$aln = $alnin->next_aln;

my $ct = 0;

# compare all to test seq

for my $ls ( sort keys %testdata ) {
    $ct++;
    my $str = $aln->bracket_string(
        -refseq  => 'testseq',
        -allele1 => 'allele1',
        -allele2 => $ls,
    );
    is( $str, $testdata{$ls}, "BIC:$str" );
}

%testdata = (
    'allele1' => 'GGATCCATT{C.C}CTACT',
    'allele2' => 'GGAT{C.-}{C.-}ATT{C.C}CT{A.C}CT',
    'allele3' => 'G{G.C}ATCCATT{C.G}CTACT',
    'allele4' => 'GGATCCATT{C.G}CTACT',
    'allele5' => 'GGATCCATT{C.G}CTAC{T.A}',
    'allele6' => 'GGATCCATT{C.G}CTA{C.G}{T.A}',
    'testseq' => 'GGATCCATT{C.G}CTACT'
);

for my $ls ( sort keys %testdata ) {
    $ct++;
    my $str = $aln->bracket_string(
        -refseq     => 'testseq',
        -allele1    => 'allele1',
        -allele2    => $ls,
        -delimiters => '{}',
        -separator  => '.'
    );
    is( $str, $testdata{$ls}, "BIC:$str" );
}

# is _remove_col really working correctly?
my $a = Bio::LocatableSeq->new(
    -id    => 'a',
    -strand => 1,
    -seq   => 'atcgatcgatcgatcg',
    -start => 5,
    -end   => 20
);
my $b = Bio::LocatableSeq->new(
    -id    => 'b',
    -strand => 1,
    -seq   => '-tcgatc-atcgatcg',
    -start => 30,
    -end   => 43
);
my $c = Bio::LocatableSeq->new(
    -id    => 'c',
    -strand => -1,
    -seq   => 'atcgatcgatc-atc-',
    -start => 50,
    -end   => 63
);
my $d = Bio::LocatableSeq->new(
    -id    => 'd',
    -strand => -1,
    -seq   => '--cgatcgatcgat--',
    -start => 80,
    -end   => 91
);
my $e = Bio::LocatableSeq->new(
    -id    => 'e',
    -strand => 1,
    -seq   => '-t-gatcgatcga-c-',
    -start => 100,
    -end   => 111
);
$aln = Bio::SimpleAlign->new();
$aln->add_seq($a);
$aln->add_seq($b);
$aln->add_seq($c);

my $gapless = $aln->remove_gaps();
foreach my $seq ( $gapless->each_seq ) {
    if ( $seq->id eq 'a' ) {
        is $seq->start, 6;
        is $seq->end,   19;
        is $seq->seq,   'tcgatcatcatc';
    }
    elsif ( $seq->id eq 'b' ) {
        is $seq->start, 30;
        is $seq->end,   42;
        is $seq->seq,   'tcgatcatcatc';
    }
    elsif ( $seq->id eq 'c' ) {
        is $seq->start, 51;
        is $seq->end,   63;
        is $seq->seq,   'tcgatcatcatc';
    }
}

$aln->add_seq($d);
$aln->add_seq($e);
$gapless = $aln->remove_gaps();
foreach my $seq ( $gapless->each_seq ) {
    if ( $seq->id eq 'a' ) {
        is $seq->start, 8;
        is $seq->end,   17;
        is $seq->seq,   'gatcatca';
    }
    elsif ( $seq->id eq 'b' ) {
        is $seq->start, 32;
        is $seq->end,   40;
        is $seq->seq,   'gatcatca';
    }
    elsif ( $seq->id eq 'c' ) {
        is $seq->start, 53;
        is $seq->end,   61;
        is $seq->seq,   'gatcatca';
    }
    elsif ( $seq->id eq 'd' ) {
        is $seq->start, 81;
        is $seq->end,   90;
        is $seq->seq,   'gatcatca';
    }
    elsif ( $seq->id eq 'e' ) {
        is $seq->start, 101;
        is $seq->end,   110;
        is $seq->seq,   'gatcatca';
    }
}

# bug 3077

my $slice = $aln->slice(1,4);

for my $seq ($slice->each_seq) {
    if ( $seq->id eq 'a' ) {
        is $seq->start, 5;
        is $seq->end,   8;
        is $seq->strand, 1;
        is $seq->seq,   'atcg';
    }
    elsif ( $seq->id eq 'b' ) {
        is $seq->start, 30;
        is $seq->end,   32;
        is $seq->strand, 1;
        is $seq->seq,   '-tcg';
    }
    elsif ( $seq->id eq 'c' ) {
        is $seq->start, 60;
        is $seq->end,   63;
        is $seq->strand, -1;
        is $seq->seq,   'atcg';
    }
    elsif ( $seq->id eq 'd' ) {
        is $seq->start, 90;
        is $seq->end,   91;
        is $seq->strand, -1;
        is $seq->seq,   '--cg';
    }
    elsif ( $seq->id eq 'e' ) {
        is $seq->start, 100;
        is $seq->end,   101;
        is $seq->strand, 1;
        is $seq->seq,   '-t-g';
    }
}

my $f = Bio::LocatableSeq->new(
    -id    => 'f',
    -seq   => 'a-cgatcgatcgat-g',
    -start => 30,
    -end   => 43
);
$aln = Bio::SimpleAlign->new();
$aln->add_seq($a);
$aln->add_seq($f);

$gapless = $aln->remove_gaps();
foreach my $seq ( $gapless->each_seq ) {
    if ( $seq->id eq 'a' ) {
        is $seq->start, 5;
        is $seq->end,   20;
        is $seq->seq,   'acgatcgatcgatg';
    }
    elsif ( $seq->id eq 'f' ) {
        is $seq->start, 30;
        is $seq->end,   43;
        is $seq->seq,   'acgatcgatcgatg';
    }
}

my $g =
  Bio::LocatableSeq->new( -id => 'g', -seq => 'atgc', -start => 5, -end => 8 );
my $h = Bio::LocatableSeq->new(
    -id    => 'h',
    -seq   => '-tcg',
    -start => 30,
    -end   => 32
);
$aln = Bio::SimpleAlign->new();
$aln->add_seq($g);
$aln->add_seq($h);

# test for new method in API get_seq_by_id
my $retrieved = $aln->get_seq_by_id('g');
is( defined $retrieved, 1 );
my $removed = $aln->remove_columns( [ 1, 3 ] );
foreach my $seq ( $removed->each_seq ) {
    if ( $seq->id eq 'g' ) {
        is $seq->start, 5;
        is $seq->end,   5;
        is $seq->seq,   'a';
    }
    elsif ( $seq->id eq 'h' ) {
        is $seq->start, 0;
        is $seq->end,   0;
        is $seq->seq,   '-';
    }
}

# work out mask_columns(), see bug 2842
SKIP: {
    test_skip(-tests => 6, -requires_module => 'IO::String');
    my $io = Bio::AlignIO->new( -file => test_input_file("testaln.clustalw") );
    my $aln = $io->next_aln();
    isa_ok( $aln, 'Bio::SimpleAlign' );
    my $consensus = <<EOU;
MNEGEHQIKLDELFEKLLRARKIFKNKDVLRHSWEPKDLPHRHEQIEALAQILV
PVLRGETMKIIFCGHHACELGEDRGTKGFVIDELKDVDEDRNGKVDVIEINCEH
MDTHYRVLPNIAKLFDDCTGIGVPMHGGPTDEVTAKLKQVIDMKERFVIIVLDE
IDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEE
EVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDL
LRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLL
DENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSK
GRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
EOU
    $consensus =~ s/\n//g;

    is( $aln->consensus_string, $consensus, 'consensus string looks ok' );


    my @cons_got = $aln->consensus_conservation;
    # 422 positions, mostly two of six sequences conserved, set as default
    my @cons_expect = (100 * 2/6) x 422;
    # Exceptionally columns as a mask, manually determined (1-based columns)
    $cons_expect[$_-1] = 100 * 1/6 for (5,12,41,70,82,310,390);
    $cons_expect[$_-1] = 100 * 3/6 for (27,30,32,36,47,49,61,66,69,71,77,79,
        81,91,96,97,105,114,115,117,118,121,122,129,140,146,156,159,160,162,
        183,197,217,221,229,242,247,248,261,266,282,287,295,316,323,329,335,337,344,);
    $cons_expect[$_-1] = 100 * 4/6 for (84,93,99,100,102,107,108,112,113,119,150,);
    $cons_expect[$_-1] = 100 * 5/6 for (81,110);
    # Format for string comparison
    @cons_expect = map { sprintf "%4.1f", $_ } @cons_expect;
    @cons_got = map { sprintf "%4.1f", $_ } @cons_got;
    is(length($aln->consensus_string), scalar(@cons_got),"conservation length");
    is_deeply(\@cons_got, \@cons_expect, "conservation scores");


    is( aln2str( $aln => 'pfam' ), <<EOA, 'looks like correct unmasked alignment (from clustalw)' );
P84139/1-420              MNEGEHQIKLDELFEKLLRARKIFKNKDVLRHSYTPKDLPLRHEQIETLAQILVPVLRGETPSNIFVYG-KTGTGKTVTVK-FVTEELKRISEKYNIPVDVIYINCEIVDTHYRVLANIVNYFKDETGIGVPMVGWPTDEVYAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTRPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLNVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
P814153/1-420             MNEGMHQIKLDVLFEKLLRARKIFKNKDVLRHSYTPKDLPHRHEQIETLAQILVPVLRGETPSNIFVYG-KTGTGKTVTVK-FVTEELKRISEKYNIPVDVIYINCEIVDTHYRVLANIVNYFKDETGIEVPMVGWPTDEVYAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
P851414/1-60              -------------------------------------------------------------MKIVWCGH-ACFLVEDRGTK-ILIDPYPDVDEDRIGKVDYILQTHEHMD-HYGKTPLIAKLSD----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
P841414/1-60              -------------------------------------------------------------MKIVWCGH-ACFLVEDRGTK-ILIDPYPDVDEDRIGKVDYILVTHEHMD-HYGKTPLIAKLSD----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
BAB68554/1-141            --------------------MLTEDDKQLIQHVWEKVLEHQEDFGAEALERMFIVYPSTKTYFPHFDLHHDSEQIRHHGKK-VVGALGDAVKHIDNLSATLSELSNLHCY-NLRVDPVNFKLLSHCFQVVLGAHLG--REYTPQVQVAYDKFLAAVSAVLAEKYR-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gb|443893|124775/1-331    -MRFRFGVVVPPAVAGARPELLVVGSRPELG-RWEPRGAVRLRPAGTAAGDGALALQEPGLWLGEVELA-AEEAAQDGAEPGRVDTFWYKFLKREPGGELSWEGNGPHHDRCCTYNENNLVDGVYCLPIG---HWGEATGHTNEMKHTTDFYFNIAGHQAMHYSRILPNIWLGSCPRQVEHVTIKLKHELGITAVMN-FQTEWDIVQNSSGCNRYPEPMTPDTMIKLYREEGLAYIWMP-TPDMSTEGRVQMLPQAVCLLHALLEKGHIVY-----VHCNAGVGRSTAAVCGWLQYVMGWNLRKVQYFLMAKRPAVYIDEEALARAQEDFFQKFGKVRSSVCSL------------------------------------------------------------------------------
EOA

    my $newaln = $aln->mask_columns(12,20,'?');
    is( aln2str( $newaln, 'pfam' ), <<EOA, 'looks like correct masked alignment (from clustalw)' );
P84139/1-420              MNEGEHQIKLD?????????RKIFKNKDVLRHSYTPKDLPLRHEQIETLAQILVPVLRGETPSNIFVYG-KTGTGKTVTVK-FVTEELKRISEKYNIPVDVIYINCEIVDTHYRVLANIVNYFKDETGIGVPMVGWPTDEVYAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTRPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLNVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
P814153/1-420             MNEGMHQIKLD?????????RKIFKNKDVLRHSYTPKDLPHRHEQIETLAQILVPVLRGETPSNIFVYG-KTGTGKTVTVK-FVTEELKRISEKYNIPVDVIYINCEIVDTHYRVLANIVNYFKDETGIEVPMVGWPTDEVYAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
P851414/1-60              -------------------------------------------------------------MKIVWCGH-ACFLVEDRGTK-ILIDPYPDVDEDRIGKVDYILQTHEHMD-HYGKTPLIAKLSD----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
P841414/1-60              -------------------------------------------------------------MKIVWCGH-ACFLVEDRGTK-ILIDPYPDVDEDRIGKVDYILVTHEHMD-HYGKTPLIAKLSD----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
BAB68554/1-141            --------------------MLTEDDKQLIQHVWEKVLEHQEDFGAEALERMFIVYPSTKTYFPHFDLHHDSEQIRHHGKK-VVGALGDAVKHIDNLSATLSELSNLHCY-NLRVDPVNFKLLSHCFQVVLGAHLG--REYTPQVQVAYDKFLAAVSAVLAEKYR-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gb|443893|124775/1-331    -MRFRFGVVVP?????????LLVVGSRPELG-RWEPRGAVRLRPAGTAAGDGALALQEPGLWLGEVELA-AEEAAQDGAEPGRVDTFWYKFLKREPGGELSWEGNGPHHDRCCTYNENNLVDGVYCLPIG---HWGEATGHTNEMKHTTDFYFNIAGHQAMHYSRILPNIWLGSCPRQVEHVTIKLKHELGITAVMN-FQTEWDIVQNSSGCNRYPEPMTPDTMIKLYREEGLAYIWMP-TPDMSTEGRVQMLPQAVCLLHALLEKGHIVY-----VHCNAGVGRSTAAVCGWLQYVMGWNLRKVQYFLMAKRPAVYIDEEALARAQEDFFQKFGKVRSSVCSL------------------------------------------------------------------------------
EOA

###### test with phylip

    my $phylip_str = <<EOF;
 3 37
seq1         AAAATGGGGG TGGT------ GGTACCT--- -------
seq2         -----GGCGG TGGTGNNNNG GGTTCCCTNN NNNNNNN
new          AAAATGGNGG TGGTN----N GGTNCCNTNN NNNNNNN

EOF

    my $phylip_masked = <<EOF;
 3 37
seq1         AAAATGGGGG TGGT------ GGTACCT--- -------
seq2         -----GGCGG TGGT?????? GGTTCCCTNN NNNNNNN
new          AAAATGGNGG TGGT?----? GGTNCCNTNN NNNNNNN

EOF

    my $phy_fh = IO::String->new( $phylip_str );

    my $in = Bio::AlignIO->new( -fh => $phy_fh, -format => 'phylip' );
    unified_diff;

    $aln = $in->next_aln();
    eq_or_diff( aln2str( $aln, 'phylip' ), $phylip_str );

    $newaln = $aln->mask_columns(15,20,'?');
    eq_or_diff( aln2str( $newaln,'phylip' ), $phylip_masked, 'align after looks ok' );
}

######## SUBROUTINES

sub aln2str {
    my ( $aln, $fmt ) = @_;
    my $out;
    my $out_fh = IO::String->new( $out );
    my $alignio_out = Bio::AlignIO->new(-fh => $out_fh, -format => $fmt);
    $alignio_out->write_aln( $aln );
    return $out;
}
