# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 30);
	
	use_ok('Bio::Seq::LargePrimarySeq');
	use_ok('Bio::Seq::LargeSeq');
	use_ok('Bio::Location::Simple');
	use_ok('Bio::Location::Fuzzy');
	use_ok('Bio::Location::Split');
	use_ok('Bio::SeqIO');
}

my $pseq = Bio::Seq::LargePrimarySeq->new();
ok $pseq;
$pseq->add_sequence_as_string('ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAAT');
$pseq->add_sequence_as_string('GTTTGGGGTTAAACCCCTTTGGGGGGT');

is $pseq->display_id('hello'), 'hello';

is $pseq->seq, 'ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAATGTTTGGGGTTAAACCCCTTTGGGGGGT' , "Sequence is " . $pseq->seq;

is $pseq->subseq(3,7), 'GGGGT', "Subseq is ".$pseq->subseq(3,7);
my $location = Bio::Location::Simple->new(-start => 4, -end => 8,
					 -strand => 1);
is($pseq->subseq($location), 'GGGTG');

my $splitlocation = Bio::Location::Split->new();

$splitlocation->add_sub_Location( Bio::Location::Simple->new('-start' => 1,
							    '-end'   => 15,
							    '-strand' => 1));

$splitlocation->add_sub_Location( Bio::Location::Simple->new('-start' => 21,
							    '-end'   => 27,
							    '-strand' => -1));

is( $pseq->subseq($splitlocation), 'ATGGGGTGGGGTGAACCCCCAA');

my $fuzzy = Bio::Location::Fuzzy->new(-start => '<10',
				     -end   => '18',
				     -strand => 1);

is( $pseq->subseq($fuzzy), 'GGTGAAACC');


is($pseq->trunc(8,15)->seq, 'GGGGTGAA',
    'trunc seq was ' . $pseq->trunc(8,15)->seq);


is $pseq->alphabet('dna'), 'dna'; # so translate will not complain
is $pseq->translate()->seq, 'MGWGETLWGWGKCLGLNPFGG';


my $seq = Bio::Seq::LargeSeq->new(-primaryseq => $pseq );

is $seq->display_id('hello'), 'hello';

is $seq->seq, 'ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAATGTTTGGGGTTAAACCCCTTTGGGGGGT' , "Sequence is " . $seq->seq;

# test SeqIO::fasta (allows LargeSeqI; bug 2490)

SKIP: {
    eval {require IO::String};
    skip "SeqIO output for LargeSeq requires IO::String", 2 if $@;
	my $str;
	my $strobj = IO::String->new($str);
	my $out = Bio::SeqIO->new(-fh => $strobj, -format => 'fasta');
    ok($out->write_seq($seq));
	like($str, qr/>hello\nATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAATGTTTGGGGTTAAACCCCTTTGGG\nGGGT/,
		 'output via Bio::SeqIO::fasta');
}						  

is $seq->subseq(3,7), 'GGGGT', "Subseq is ".$seq->subseq(3,7);
is ($seq->trunc(8,15)->seq, 'GGGGTGAA', 
    'trunc seq was ' . $seq->trunc(8,15)->seq);

is $seq->alphabet('dna'), 'dna'; # so translate will not complain
is $seq->translate()->seq, 'MGWGETLWGWGKCLGLNPFGG';

$seq = Bio::Seq::LargeSeq->new( -display_id => 'hello');
$seq->seq('ATGGGGTGGGGT');
is $seq->display_id, 'hello';

is $seq->seq, 'ATGGGGTGGGGT' , "Sequence is " . $seq->seq;

is $seq->subseq(3,7), 'GGGGT', "Subseq is ".$seq->subseq(3,7);
is ($seq->trunc(8,12)->seq, 'GGGGT', 
    'trunc seq was ' . $seq->trunc(8,12)->seq);

is $seq->alphabet('dna'), 'dna'; # so translate will not complain
is $seq->translate()->seq, 'MGWG';
