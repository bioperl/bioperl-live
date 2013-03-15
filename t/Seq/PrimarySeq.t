# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use Data::Dumper;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin( -tests => 185 );

    use_ok('Bio::PrimarySeq');
    use_ok('Bio::Location::Simple');
    use_ok('Bio::Location::Fuzzy');
    use_ok('Bio::Location::Split');
}


# Bare object
ok my $seq = Bio::PrimarySeq->new(), 'Bare object';
isa_ok $seq, 'Bio::PrimarySeqI';
is $seq->id, undef;
is $seq->seq, undef;
is $seq->length, 0;
is $seq->alphabet, undef;
is $seq->is_circular, undef;


# Empty sequence
ok $seq = Bio::PrimarySeq->new( -seq => '', -nowarnonempty => 1);
is $seq->seq, '';
is $seq->length, 0;
is $seq->alphabet, undef;


# Basic tests
ok $seq = Bio::PrimarySeq->new(
    '-seq'              => 'TTGGTGGCGTCAACT',
    '-display_id'       => 'new-id',
    '-alphabet'         => 'dna',
    '-accession_number' => 'X677667',
    '-desc'             => 'Sample Bio::Seq object'
);
ok defined $seq;
is $seq->accession_number(), 'X677667';
is $seq->seq(),              'TTGGTGGCGTCAACT';
is $seq->display_id(),       'new-id';
is $seq->alphabet(),         'dna';
is $seq->is_circular(),      undef;
ok $seq->is_circular(1);
is $seq->is_circular(0), 0;

# check IdentifiableI and DescribableI interfaces
isa_ok $seq, 'Bio::IdentifiableI';
isa_ok $seq, 'Bio::DescribableI';

# make sure all methods are implemented
is $seq->authority("bioperl.org"), "bioperl.org";
is $seq->authority, "bioperl.org";
is $seq->namespace("t"), "t";
is $seq->namespace, "t";
is $seq->version(0), 0;
is $seq->version, 0;
is $seq->lsid_string(), "bioperl.org:t:X677667";
is $seq->namespace_string, "t:X677667.0";
is $seq->version(47), 47;
is $seq->version, 47;
is $seq->namespace_string, "t:X677667.47";
is $seq->description, 'Sample Bio::Seq object';
is $seq->display_name, "new-id";


# Test subseq
is $seq->subseq(2, 5), 'TGGT';

is $seq->subseq( -start => 1, -end => 15), 'TTGGTGGCGTCAACT';

my $location = Bio::Location::Simple->new(
    '-start'  => 2,
    '-end'    => 5,
    '-strand' => -1
);
is $seq->subseq($location), 'ACCA';

my $splitlocation = Bio::Location::Split->new();
$splitlocation->add_sub_Location(
    Bio::Location::Simple->new(
        '-start'  => 1,
        '-end'    => 4,
        '-strand' => 1
    )
);

$splitlocation->add_sub_Location(
    Bio::Location::Simple->new(
        '-start'  => 7,
        '-end'    => 12,
        '-strand' => -1
    )
);

is $seq->subseq($splitlocation), 'TTGGTGACGC';

my $fuzzy = Bio::Location::Fuzzy->new(
    -start  => '<3',
    -end    => '8',
    -strand => 1
);

is $seq->subseq($fuzzy), 'GGTGGC';

{
    ok my $seq = Bio::PrimarySeq->new( -seq => 'TT-GTGGCGTCAACT' );
    is $seq->subseq(2, 5, 'nogap'), 'TGT';
    is $seq->subseq( -start => 2, -end => 5, -nogap => 1 ), 'TGT';
    my $location = Bio::Location::Simple->new(
       '-start'  => 2,
       '-end'    => 5,
       '-strand' => 1
    );
    is $seq->subseq( $location, -nogap => 1), 'TGT';

    is $seq->subseq(-start=>2, -end=>5, -replace_with=>'aa'), 'T-GT';
    is $seq->seq, 'TaaGGCGTCAACT';

    throws_ok { $seq->subseq(-start=>2, -end=>5, -replace_with=>'?!'); } qr/.+/;
}

{
    ok my $seq = Bio::PrimarySeq->new( -seq => 'AACCGGTT', -is_circular => 1 );
    is $seq->subseq( -start => 7, -end => 10 ), 'TTAA';
}

### Test for Bug #2936
# Without strand input argument (case: user don't think is necessary)
my $split_loc_obj1 = Bio::Location::Split->new();
$split_loc_obj1->add_sub_Location(
    Bio::Location::Simple->new(
        '-start'  => 1,
        '-end'    => 10
    )
);
$split_loc_obj1->add_sub_Location(
    Bio::Location::Simple->new(
        '-start'  => 20,
        '-end'    => 30
    )
);
# With strand input argument (case: user provides the argument)
my $split_loc_obj2 = Bio::Location::Split->new();
$split_loc_obj2->add_sub_Location(
    Bio::Location::Simple->new(
        '-start'  => 1,
        '-end'    => 10,
        '-strand' => 1
    )
);
$split_loc_obj2->add_sub_Location(
    Bio::Location::Simple->new(
        '-start'  => 20,
        '-end'    => 30,
        '-strand' => 1
    )
);
is $split_loc_obj1->to_FTstring, "join(1..10,20..30)";
is $split_loc_obj2->to_FTstring, "join(1..10,20..30)";
$split_loc_obj1->flip_strand;
$split_loc_obj2->flip_strand;
is $split_loc_obj1->to_FTstring, "complement(join(1..10,20..30))";
is $split_loc_obj2->to_FTstring, "complement(join(1..10,20..30))";
###

# Test trunc
my $trunc = $seq->trunc( 1, 4 );
isa_ok $trunc, 'Bio::PrimarySeqI';
is $trunc->seq(), 'TTGG' or diag( "Expecting TTGG. Got " . $trunc->seq() );

$trunc = $seq->trunc($splitlocation);
isa_ok $trunc, 'Bio::PrimarySeqI' ;
is $trunc->seq(), 'TTGGTGACGC';

$trunc = $seq->trunc($fuzzy);
isa_ok $trunc, 'Bio::PrimarySeqI';
is $trunc->seq(), 'GGTGGC';

my $rev = $seq->revcom();
isa_ok $rev, 'Bio::PrimarySeqI';

is $rev->seq(), 'AGTTGACGCCACCAA'
  or diag( 'revcom() failed, was ' . $rev->seq() );

is $rev->display_id, 'new-id';
is $rev->display_name(), 'new-id';
is $rev->accession_number(), 'X677667';
is $rev->alphabet, 'dna';
is $rev->description, 'Sample Bio::Seq object';


TODO: {
    local $TODO =
      'all attributes of primaryseqs are not currently copied through revcom()';
    # Probably also not copied through trunc(), transcribe() and rev_transcribe()
    is $rev->is_circular(), 0,           'is_circular copied through revcom';
    is $rev->version, 47,                'version copied through revcom';
    is $rev->authority, 'bioperl.org',   'authority copied through revcom';
    is $rev->namespace, 't',             'namespace copied through revcom';
    is $rev->namespace_string(),
        "t:X677667.47", 'namespace_string copied through revcom';
}

#
# Translate
#

my $aa = $seq->translate();    # TTG GTG GCG TCA ACT
is $aa->seq, 'LVAST', "Translation: " . $aa->seq;

# tests for non-standard initiator codon coding for
# M by making translate() look for an initiator codon and
# terminator codon ("complete", the 5th argument below)
$seq->seq('TTGGTGGCGTCAACTTAA');    # TTG GTG GCG TCA ACT TAA
$aa = $seq->translate( undef, undef, undef, undef, 1 );
is $aa->seq, 'MVAST', "Translation: " . $aa->seq;

# same test as previous, but using named parameter
$aa = $seq->translate( -complete => 1 );
is $aa->seq, 'MVAST', "Translation: " . $aa->seq;

# find ORF, ignore codons outside the ORF or CDS
$seq->seq('TTTTATGGTGGCGTCAACTTAATTT');    # ATG GTG GCG TCA ACT
$aa = $seq->translate( -orf => 1 );
is $aa->seq, 'MVAST*', "Translation: " . $aa->seq;

# smallest possible ORF
$seq->seq("ggggggatgtagcccc");             # atg tga
$aa = $seq->translate( -orf => 1 );
is $aa->seq, 'M*', "Translation: " . $aa->seq;

# same as previous but complete, so * is removed
$aa = $seq->translate(
    -orf      => 1,
    -complete => 1
);
is $aa->seq, 'M', "Translation: " . $aa->seq;

# ORF without termination codon
# should warn, let's change it into throw for testing
$seq->verbose(2);
$seq->seq("ggggggatgtggcccc");    # atg tgg ccc
eval { $seq->translate( -orf => 1 ); };
like( $@, qr/\batgtggccc\b/i );
$seq->verbose(-1);
$aa = $seq->translate( -orf => 1 );
is $aa->seq, 'MWP', "Translation: MWP";
$seq->verbose(0);

# use non-standard codon table where terminator is read as Q
$seq->seq('ATGGTGGCGTCAACTTAG');    # ATG GTG GCG TCA ACT TAG
$aa = $seq->translate( -codontable_id => 6 );
is $aa->seq, 'MVASTQ' or diag( "Translation: " . $aa->seq );

# insert an odd character instead of terminating with *
$aa = $seq->translate( -terminator => 'X' );
is $aa->seq, 'MVASTX' or diag( "Translation: " . $aa->seq );

# change frame from default
$aa = $seq->translate( -frame => 1 );    # TGG TGG CGT CAA CTT AG
is $aa->seq, 'WWRQL' or diag( "Translation: " . $aa->seq );

$aa = $seq->translate( -frame => 2 );    # GGT GGC GTC AAC TTA G
is $aa->seq, 'GGVNL' or diag( "Translation: " . $aa->seq );

# TTG is initiator in Standard codon table? Afraid so.
$seq->seq("ggggggttgtagcccc");           # ttg tag
$aa = $seq->translate( -orf => 1 );
is $aa->seq, 'L*' or diag( "Translation: " . $aa->seq );

# Replace L at 1st position with M by setting complete to 1
$seq->seq("ggggggttgtagcccc");           # ttg tag
$aa = $seq->translate(
    -orf      => 1,
    -complete => 1
);
is $aa->seq, 'M' or diag( "Translation: " . $aa->seq );

# Ignore non-ATG initiators (e.g. TTG) in codon table
$seq->seq("ggggggttgatgtagcccc");        # atg tag
$aa = $seq->translate(
    -orf      => 1,
    -start    => "atg",
    -complete => 1
);
is $aa->seq, 'M' or diag( "Translation: " . $aa->seq );

# test for character '?' in the sequence string
is $seq->seq('TTGGTGGCG?CAACT'), 'TTGGTGGCG?CAACT';

# test for some aliases
$seq = Bio::PrimarySeq->new(
    -id          => 'aliasid',
    -description => 'Alias desc'
);
is $seq->description, 'Alias desc';
is $seq->display_id,  'aliasid';

# Test alphabet

ok $seq->seq('actgx');
is $seq->alphabet, 'protein', 'Alphabet';
ok $seq->seq('actge');
is $seq->alphabet, 'protein';
ok $seq->seq('actgf');
is $seq->alphabet, 'protein';
ok $seq->seq('actgi');
is $seq->alphabet, 'protein';
ok $seq->seq('actgj');
is $seq->alphabet, 'protein';
ok $seq->seq('actgl');
is $seq->alphabet, 'protein';
ok $seq->seq('actgo');
is $seq->alphabet, 'protein';
ok $seq->seq('actgp');
is $seq->alphabet, 'protein';
ok $seq->seq('actgq');
is $seq->alphabet, 'protein';
ok $seq->seq('actgz');
is $seq->alphabet, 'protein';
ok $seq->seq('actgn');
is $seq->alphabet, 'dna';
ok $seq->seq('acugn');
is $seq->alphabet, 'rna';
ok $seq->seq('bdhkm');
is $seq->alphabet, 'protein';
ok $seq->seq('rsvwx');
is $seq->alphabet, 'protein';
ok $seq->seq('AAACTYAAAAGAATTGRCGG'); # valid degenerate DNA PCR primer sequence (90% ACGTN)
is $seq->alphabet, 'dna';
ok $seq->seq('AAACTYAAAKGAATTGRCGG'); # another primer previously detected as protein (85% ACGTN)
is $seq->alphabet, 'dna';
ok $seq->seq('YWACTYAAAKGARTTGRCGG'); # 70% ACGTNWSRM. Everything <= 70% is considered a protein
is $seq->alphabet, 'dna';
ok $seq->seq('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'); # Bug 2438
is $seq->alphabet, 'protein', 'Bug 2438';
ok $seq->seq('CAGTCXXXXXXXXXXXXXXXXXXXXXXXXXXXCAGCG');
is $seq->alphabet, 'protein';
ok $seq->seq('WTGGGGCTATGAAAAAAAAAWTTKMGMMAAAAAWTTWTKRWMRATC'); # showed up on MAKER list
is $seq->alphabet, 'dna';

ok $seq->seq('actgn', 'protein'); # accept specified alphabet, no matter what
is $seq->alphabet, 'protein';
ok $seq->seq('bdhkm', 'dna');
is $seq->alphabet, 'dna';


# Bug #2864:

$seq = Bio::PrimarySeq->new( -display_id => 0, -seq => 'GATC' );

is $seq->display_id, 0, "Bug #2864";

# Test that the check for terminators inside the translated protein
# works when the terminator isn't '*':

$seq = Bio::PrimarySeq->new(-seq=>'ATGCTCTAAGCAGGGTAA'); # ML*AG*
eval { $aa = $seq->translate(-complete=>1, -throw=>1, -terminator=>'#') };
my $error = $@;
ok $error =~ /\QTerminator codon inside CDS!\E/, 'Terminator + inside sequence';

$seq = Bio::PrimarySeq->new(-seq=>'ATGCTCGCAGGGTAA'); # MLAG*
$aa = $seq->translate(-complete=>1, -throw=>1, -terminator=>'#');
is $aa->seq, 'MLAG';


# Test length method
ok $seq = Bio::PrimarySeq->new(), 'Length method';
is $seq->length, 0;
ok $seq->length(123);
is $seq->length, 123;

ok $seq = Bio::PrimarySeq->new( -seq => 'ATGCTCTAAGCAGGGTAA' );
is $seq->length, 18;
ok $seq->seq('ATGCTCTAAG');
is $seq->length, 10;
is $seq->seq(undef), undef;
is $seq->length, 0;

ok $seq = Bio::PrimarySeq->new( -length => 123 );
is $seq->length, 123;

ok $seq = Bio::PrimarySeq->new( -seq => 'ATGCTCTAAGCAGGGTAA' );
is $seq->length, 18;
ok $seq->length( $seq->length ); # save memory by removing seq
is $seq->seq( undef ), undef;    # ... but keeping a record of length
is $seq->length, 18;
is $seq->seq, undef;
ok $seq->seq('ACGT');
is $seq->length, 4; # manually-specified length changed when sequence is changed

throws_ok { $seq->length(666); } qr/.+/; # Cannot lie about length


# Sequence validation method
is $seq->validate_seq( undef    ), 1;
is $seq->validate_seq( ''       ), 1;
is $seq->validate_seq( 'acgt'   ), 1;
is $seq->validate_seq( 'ACGT'   ), 1;
is $seq->validate_seq( 'XFRH'   ), 1;
is $seq->validate_seq( '-~'     ), 1; # gap symbols
is $seq->validate_seq( '-.*?=~' ), 1; # other valid symbols
is $seq->validate_seq( '0'      ), 0;
is $seq->validate_seq( '   '    ), 0;
is $seq->validate_seq( 'AAAA$'  ), 0;
is $seq->validate_seq( 'tt&t!'  ), 0;

throws_ok { $seq->validate_seq('tt&t!', 1); } qr/.+/;


# Test direct option (no sequence validation)
throws_ok { $seq = Bio::PrimarySeq->new(-seq => 'A\T$AGQ+T'); } qr/.+/, 'Validation';
ok $seq = Bio::PrimarySeq->new( -seq => 'A\T$AGQ+T', -direct => 1 );
is $seq->seq, 'A\T$AGQ+T';
throws_ok { $seq->seq('NT@/') } qr/.+/;

# Set a sequence by reference
my $string = 'AAAACCCCGGGGTTTT';
ok $seq = Bio::PrimarySeq->new( -ref_to_seq => \$string );
is $seq->seq, 'AAAACCCCGGGGTTTT';


# Test internal PrimarySeqI _find_orfs function and translate( -orf => 'longest' )
{
    my @tests = (
        #tiny test
        ['TTTTATGGTGGCGTCAACTTAATTT',
         [[4,22,18,1]],
        ],

        #bigger test (this is a tomato unigene)
        ['GAAGGCTGGTTCTGAGTTGGATCTATGTTTGATGAAGGGAAGTAGACCGGAGGTCTTGCATCAGCAATATTAGTACCAAATCCAGGTGGAGGCGCATCCTGTCTCCGTTGCATTTCAACTTTCATTTCAGCAATCTGTTGCATCAGTTGCATGATCAATTCATTCTGTTCCACTACAGTGGGCTGAGCGACCACAACGTCAGTAAGACGCCCTTCGTCATTGTTGTCTCCCATAACTGTTTTTCCTTTATCTGAATTTGATCGAGGGAAGGAATCTGTAGGACCTTTCGATCTGGTGAAGTAAGGATGATCTGCCAGCTTTATTGACACAGATCAGTAAAAAGGTACCTGAAAGGTAAAAACAACTCAAAGGCAAATTTGTTAGTGCATATCCAGAGTACAAAATGCTTAATATCGCACATAAAACCGATAAACACACAAGTCGTTTTGTTTGAGGATATCTTAACCCACGAATAAGGACGGATATATATTTTGAACAAACAGGAATTTGTTTGTTTGGCGTTATCTTGGGAAATCTG',
         [[98,254,156,2],[347,476,129,2],[219,303,84,0],[16,73,57,1],[403,454,51,1],[310,358,48,1],[235,280,45,1],[491,536,45,2],[150,186,36,0],[507,537,30,0],[5,32,27,2],[511,538,27,1],[24,45,21,0],[305,326,21,2],[450,465,15,0]],
        ],


       );
    foreach my $test (@tests) {
        my ($test_seq, $orfs) = @$test;
        my @orfs = Bio::PrimarySeqI::_find_orfs_nucleotide(
            undef,
            $test_seq,
            Bio::Tools::CodonTable->new,
            undef,
           ); # ATG GTG GCG TCA ACT
        is_deeply( \@orfs, $orfs, '_find_orfs 1')
            or diag "for $test_seq, _find_orfs returned:\n"
                    .Dumper([map [@$_], @orfs]);

        is_deeply( $orfs->[0],
                   (sort {$b->[2] <=> $a->[2]} @$orfs)[0],
                   'orfs are sorted by descending length'
                  );

        # make sure we get the same sequence by taking the longest orf
        # nucleotide from the test data and translating it, as by
        # calling translate with -orf => 'longest'
        is(
            Bio::PrimarySeq
              ->new( -seq => $test_seq, -id => 'fake_id' )
              ->translate( -orf => 'longest' )
              ->seq,

            Bio::PrimarySeq
              ->new( -seq => substr( $test_seq, $orfs->[0][0], $orfs->[0][2] ),
                     -id => 'foo'
                    )
              ->translate
              ->seq,
            'got correct -orf => "longest" seq',
           );
    }
}
