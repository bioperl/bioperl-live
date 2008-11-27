# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests               => 14,
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
    open(FILE, "< $filename") or die("cannot open $filename");
    my @datain = <FILE>;
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
