# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests               => 22,
               -requires_modules    => [],
               -requires_networking => 0,
              );

    use_ok('Bio::SeqIO::fasta');
}

my $verbose = test_debug();

my $format = 'fasta';
my $seqio_obj = Bio::SeqIO->new(-file   => test_input_file("test.$format"),
                                -format => $format);

isa_ok($seqio_obj, 'Bio::SeqIO');

my @methods = qw(next_seq write_seq);
foreach my $method (@methods) {
    can_ok($seqio_obj, $method) ||
        diag "$method method not implemented for $format";
}

# checking the first sequence object
my $seq_obj = $seqio_obj->next_seq();
isa_ok($seq_obj, 'Bio::Seq');
my %expected = ('seq'         => 'MVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGL' .
                                 'DYRTTDENLKAHEKWGNIVDVVVMKDPRTKRSRGFGFI' .
                                 'TYSHSSMIDEAQKSRPHKIDGRVEPKRAVPRQDIDSPN' .
                                 'AGATVKKLFVGALKDDHDEQSIRDYFQHFGNIVDNIVI' .
                                 'DKETGKKRGFAFVEFDDYDPVDKVVLQKQHQLNGKMVD' .
                                 'VKKALPKNDQQGGGGGRGGPGGRAGGNRGNMGGGNYGN' .
                                 'QNGGGNWNNGGNNWGNNRGNDNWGNNSFGGGGGGGGGY' .
                                 'GGGNNSWGNNNPWDNGNGGGNFGGGGNNWNGGNDFGGY' .
                                 'QQNYGGGPQRGGGNFNNNRMQPYQGGGGFKAGGGNQGN' .
                                 'YGNNQGFNNGGNNRRY',
                'length'      => '358',
                'primary_id'  => 'roa1_drome',
                'description' => qr(Rea guano receptor type III),
               );
is   ($seq_obj->seq(),         $expected{'seq'},         'sequence');
is   ($seq_obj->length(),      $expected{'length'},      'length');
is   ($seq_obj->primary_id(),  $expected{'primary_id'},  'primary_id');
like ($seq_obj->description(), $expected{'description'}, 'description');


# checking the second sequence object
my $seq_obj2  = $seqio_obj->next_seq();
isa_ok($seq_obj2, 'Bio::Seq');
my %expected2 = ('seq'         => 'MVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGL' .
                                  'DYRTTDENLKAHEKWGNIVDVVVMKDPTSTSTSTSTST' .
                                  'STSTSTMIDEAQKSRPHKIDGRVEPKRAVPRQDIDSPN' .
                                  'AGATVKKLFVGALKDDHDEQSIRDYFQHLLLLLLLDLL' .
                                  'LLDLLLLDLLLFVEFDDYDPVDKVVLQKQHQLNGKMVD' .
                                  'VKKALPKNDQQGGGGGRGGPGGRAGGNRGNMGGGNYGN' .
                                  'QNGGGNWNNGGNNWGNNRGNDNWGNNSFGGGGGGGGGY' .
                                  'GGGNNSWGNNNPWDNGNGGGNFGGGGNNWNGGNDFGGY' .
                                  'QQNYGGGPQRGGGNFNNNRMQPYQGGGGFKAGGGNQGN' .
                                  'YGNNQGFNNGGNNRRY',
                 'length'      => '358',
                 'primary_id'  => 'roa2_drome',
                 'description' => qr(Rea guano ligand),
                );
is   ($seq_obj2->seq(),         $expected2{'seq'},         'sequence');
is   ($seq_obj2->length(),      $expected2{'length'},      'length');
is   ($seq_obj2->primary_id(),  $expected2{'primary_id'},  'primary_id');
like ($seq_obj2->description(), $expected2{'description'}, 'description');

# from testformats.pl
SKIP: {
    test_skip(-tests => 4, -requires_modules => [qw(Algorithm::Diff
                                                    IO::ScalarArray
                                                    IO::String)]);
    use_ok('Algorithm::Diff');
    eval "use Algorithm::Diff qw(diff LCS);";
    use_ok('IO::ScalarArray');
    use_ok('IO::String');
    
    my ($file, $type) = ("test.$format", $format);
    my $filename = test_input_file($file);
    print "processing file $filename\n" if $verbose;
    open my $FILE, '<', $filename or die "Could not read file '$filename': $!\n";
    my @datain = <$FILE>;
    close $FILE;

    my $in = new IO::String(join('', @datain));
    my $seqin = new Bio::SeqIO( -fh => $in,
                -format => $type);
    my $out = new IO::String;
    my $seqout = new Bio::SeqIO( -fh => $out,
                 -format => $type);
    my $seq;
    while( defined($seq = $seqin->next_seq) ) {
    $seqout->write_seq($seq);
    }
    $seqout->close();
    $seqin->close();
    my $strref = $out->string_ref;
    my @dataout = map { $_."\n"} split(/\n/, $$strref );
    my @diffs = &diff( \@datain, \@dataout);
    is(@diffs, 0, "$format format can round-trip");
    
    if(@diffs && $verbose) {
        foreach my $d ( @diffs ) {
            foreach my $diff ( @$d ) {
                chomp($diff->[2]);
                print $diff->[0], $diff->[1], "\n>", $diff->[2], "\n";
            }
        }
        print "in is \n", join('', @datain), "\n";
        print "out is \n", join('',@dataout), "\n";
    }

}

# bug 1508
# test genbank, gcg, ace against fasta (should throw an exception on each)

for my $file (qw(roa1.genbank test.gcg test.ace test.raw)) {
    my $in = Bio::SeqIO->new(-file   => test_input_file($file),
                             -format => 'fasta');
    throws_ok {$in->next_seq}
        qr/The sequence does not appear to be FASTA format/, "dies with $file";
}
