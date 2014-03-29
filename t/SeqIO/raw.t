# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests               => 25,
			   -requires_modules    => [],
			   -requires_networking => 0,
			  );
	
	use_ok('Bio::SeqIO::raw');
}

my $verbose = test_debug();

my $format = 'raw';
my $seqio_obj = Bio::SeqIO->new(-file   => test_input_file("test.$format"),
						        -format => $format);

isa_ok($seqio_obj, 'Bio::SeqIO');


is $seqio_obj->variant, 'multiple';

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
			   );
is   ($seq_obj->seq(),         $expected{'seq'},         'sequence');
is   ($seq_obj->length(),      $expected{'length'},      'length');


# checking the second sequence object
my $seq_obj2  = $seqio_obj->next_seq();
isa_ok($seq_obj2, 'Bio::Seq');
my %expected2 = ('seq'         => 'MNKQRGTYSEVSLAQDPKRQQRKLKGNKISISGTKQEI' .
								  'FQVELNLQNASSDHQGNDKTYHCKGLLPPPEKLTAEVL' .
								  'GIICIVLMATVLKTIVLIPCIGVLEQNNFSLNRRMQKA' .
								  'RHCGHCPEEWITYSNSCYYIGKERRTWEERVCWPVLRR' .
								  'TLICFL',
				 'length'      => '158',
			    );
is   ($seq_obj2->seq(),         $expected2{'seq'},         'sequence');
is   ($seq_obj2->length(),      $expected2{'length'},      'length');
	
# from testformats.pl
SKIP: {

    test_skip(-tests => 2,
              -requires_modules => [qw(Algorithm::Diff
                                    IO::ScalarArray
                                    IO::String)]);
    use_ok('Algorithm::Diff', qw(diff LCS));
    
	my ($file, $type) = ("test.$format", $format);
    my $filename = test_input_file($file);
    print "processing file $filename\n" if $verbose;
    open my $FILE, '<', $filename or die "Could not read file '$filename': $!\n";
    my @datain = <$FILE>;
    close $FILE;

    my $in = IO::String->new(join('', @datain));
    my $seqin = Bio::SeqIO->new( -fh => $in,
                -format => $type);
    my $out = IO::String->new;
    my $seqout = Bio::SeqIO->new( -fh => $out,
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
        use Data::Dumper;  # should be safe for 5.005 and greater
        foreach my $d ( @diffs ) {
            print STDERR Dumper $d;
            foreach my $diff ( @$d ) {
                chomp($diff->[2]);
                print $diff->[0], $diff->[1], "\n>", $diff->[2], "\n";
            }
        }
        print "in is \n", join('', @datain), "\n";
        print "out is \n", join('',@dataout), "\n";	
    }
}

# test raw variants
my @seq = qw(MVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGLDYRTTDENLKAHEKWGNIVDVV
VMKDPRTKRSRGFGFITYSHSSMIDEAQKSRPHKIDGRVEPKRAVPRQDIDSPNAGATVK
KLFVGALKDDHDEQSIRDYFQHFGNIVDNIVIDKETGKKRGFAFVEFDDYDPVDKVVLQK
QHQLNGKMVDVKKALPKNDQQGGGGGRGGPGGRAGGNRGNMGGGNYGNQNGGGNWNNGGN
NWGNNRGNDNWGNNSFGGGGGGGGGYGGGNNSWGNNNPWDNGNGGGNFGGGGNNWNGGND
FGGYQQNYGGGPQRGGGNFNNNRMQPYQGGGGFKAGGGNQGNYGNNQGFNNGGNNRRY);

$seqio_obj = Bio::SeqIO->new(-file   => test_input_file("test2.raw"),
						        -format => 'raw');

is($seqio_obj->variant, 'multiple');

my $ct = 0;

while (my $seq = $seqio_obj->next_seq) {
    is($seq->seq, $seq[$ct]);
    $ct++;
}

is($ct, 6);

$seqio_obj = Bio::SeqIO->new(-file   => test_input_file("test2.raw"),
						        -format => 'raw-single');

is($seqio_obj->variant, 'single');

my $seq = $seqio_obj->next_seq;
is($seq->seq, join('', @seq));

$seqio_obj = Bio::SeqIO->new(-file      => test_input_file("test2.raw"),
                            -format     => 'raw',
                            -variant    => 'single');

is($seqio_obj->variant, 'single');

$seq = $seqio_obj->next_seq;
is($seq->seq, join('', @seq));
