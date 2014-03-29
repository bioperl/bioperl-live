# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests               => 17,
			   -requires_modules    => [],
			   -requires_networking => 0,
			  );
	
	use_ok('Bio::SeqIO::gcg');
}

my $verbose = test_debug();

my $format = 'gcg';
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
isa_ok($seq_obj, 'Bio::Seq::RichSeq');
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
TODO: {
	local $TODO = 'possible bug: RichSeq not setting primary_id?';
	is   ($seq_obj->primary_id(),  $expected{'primary_id'},  'primary_id');
}
like ($seq_obj->description(), $expected{'description'}, 'description');

# test DOS linefeeds in gcg parser
my $str = Bio::SeqIO->new(-file => test_input_file('test_badlf.gcg'),
								  -verbose => $verbose,
								  -format => 'GCG');
ok($str);
my $seq = $str->next_seq();
isa_ok ($seq, 'Bio::SeqI');
is(length($seq->seq), $seq->length);
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n" if( $verbose);



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
