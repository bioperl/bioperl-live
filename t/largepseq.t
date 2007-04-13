# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {     
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 27;
	use_ok('Bio::Seq::LargePrimarySeq');
	use_ok('Bio::Seq::LargeSeq');
	use_ok('Bio::Location::Simple');
	use_ok('Bio::Location::Fuzzy');
	use_ok('Bio::Location::Split');
}

my $pseq = Bio::Seq::LargePrimarySeq->new();
ok $pseq;
$pseq->add_sequence_as_string('ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAAT');
$pseq->add_sequence_as_string('GTTTGGGGTTAAACCCCTTTGGGGGGT');

is $pseq->display_id('hello'), 'hello';

is $pseq->seq, 'ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAATGTTTGGGGTTAAACCCCTTTGGGGGGT' , "Sequence is " . $pseq->seq;

is $pseq->subseq(3,7), 'GGGGT', "Subseq is ".$pseq->subseq(3,7);
my $location = new Bio::Location::Simple(-start => 4, -end => 8,
					 -strand => 1);
is($pseq->subseq($location), 'GGGTG');

my $splitlocation = new Bio::Location::Split;

$splitlocation->add_sub_Location( new Bio::Location::Simple('-start' => 1,
							    '-end'   => 15,
							    '-strand' => 1));

$splitlocation->add_sub_Location( new Bio::Location::Simple('-start' => 21,
							    '-end'   => 27,
							    '-strand' => -1));

is( $pseq->subseq($splitlocation), 'ATGGGGTGGGGTGAACCCCCAA');

my $fuzzy = new Bio::Location::Fuzzy(-start => '<10',
				     -end   => '18',
				     -strand => 1);

is( $pseq->subseq($fuzzy), 'GGTGAAACC');


is($pseq->trunc(8,15)->seq, 'GGGGTGAA', 
    'trunc seq was ' . $pseq->trunc(8,15)->seq);


is $pseq->alphabet('dna'), 'dna'; # so translate will not complain
is $pseq->translate()->seq, 'MGWGETLWGWGKCLGLNPFGG';


my $seq = new Bio::Seq::LargeSeq(-primaryseq => $pseq );

is $seq->display_id('hello'), 'hello';

is $seq->seq, 'ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAATGTTTGGGGTTAAACCCCTTTGGGGGGT' , "Sequence is " . $seq->seq;

is $seq->subseq(3,7), 'GGGGT', "Subseq is ".$seq->subseq(3,7);
is ($seq->trunc(8,15)->seq, 'GGGGTGAA', 
    'trunc seq was ' . $seq->trunc(8,15)->seq);

is $seq->alphabet('dna'), 'dna'; # so translate will not complain
is $seq->translate()->seq, 'MGWGETLWGWGKCLGLNPFGG';

$seq = new Bio::Seq::LargeSeq( -display_id => 'hello');
$seq->seq('ATGGGGTGGGGT');
is $seq->display_id, 'hello';

is $seq->seq, 'ATGGGGTGGGGT' , "Sequence is " . $seq->seq;

is $seq->subseq(3,7), 'GGGGT', "Subseq is ".$seq->subseq(3,7);
is ($seq->trunc(8,12)->seq, 'GGGGT', 
    'trunc seq was ' . $seq->trunc(8,12)->seq);

is $seq->alphabet('dna'), 'dna'; # so translate will not complain
is $seq->translate()->seq, 'MGWG';
