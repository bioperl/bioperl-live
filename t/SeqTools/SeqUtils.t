# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
#    use List::MoreUtils qw(uniq);
    use Bio::Root::Test;
    
    test_begin(-tests => 53);
	
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::SeqUtils');
	use_ok('Bio::LiveSeq::Mutation');
	use_ok('Bio::SeqFeature::Generic');
	use_ok('Bio::Annotation::SimpleValue');
}

my ($seq, $util, $ascii, $ascii_aa, $ascii3);

# Entire alphabet now IUPAC-endorsed and used in GenBank (Oct 2006)          
$ascii =    'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
$ascii_aa = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

$ascii3 = 
    'AlaAsxCysAspGluPheGlyHisIleXleLysLeuMetAsnPylProGlnArgSerThrSecValTrpXaaTyrGlx';

$seq = Bio::PrimarySeq->new('-seq'=> $ascii,
			    '-alphabet'=>'protein', 
			       '-id'=>'test');

# one letter amino acid code to three letter code
ok $util = Bio::SeqUtils->new();
is $util->seq3($seq), $ascii3;

#using anonymous hash
is (Bio::SeqUtils->seq3($seq), $ascii3); 
is (Bio::SeqUtils->seq3($seq, undef, ','), 
    'Ala,Asx,Cys,Asp,Glu,Phe,Gly,His,Ile,Xle,Lys,'.
    'Leu,Met,Asn,Pyl,Pro,Gln,Arg,Ser,Thr,Sec,Val,Trp,Xaa,Tyr,Glx');

$seq->seq('asd-KJJK-');
is (Bio::SeqUtils->seq3($seq, '-', ':'), 
    'Ala:Ser:Asp:Ter:Lys:Xle:Xle:Lys:Ter');

# three letter amino acid code to one letter code
ok (Bio::SeqUtils->seq3in($seq, 'AlaPYHCysAspGlu')); 
is $seq->seq, 'AXCDE';
is (Bio::SeqUtils->seq3in($seq, $ascii3)->seq, $ascii_aa);

#
# Tests for multiframe translations
#

$seq = Bio::PrimarySeq->new('-seq'=> 'agctgctgatcggattgtgatggctggatggcttgggatgctgg',
			    '-alphabet'=>'dna', 
			    '-id'=>'test2');

my @a = $util->translate_3frames($seq);
is scalar @a, 3;
#foreach $a (@a) {
#    print 'ID: ', $a->id, ' ', $a->seq, "\n";
#}

@a = $util->translate_6frames($seq);
is scalar @a, 6;
#foreach $a (@a) {
#    print 'ID: ', $a->id, ' ', $a->seq, "\n";
#}

#
# test for valid AA return
#

my @valid_aa = sort Bio::SeqUtils->valid_aa;
is(@valid_aa, 27);
is($valid_aa[1], 'A');

@valid_aa = sort Bio::SeqUtils->valid_aa(1);
is(@valid_aa, 27);
is ($valid_aa[1], 'Arg');

my %valid_aa = Bio::SeqUtils->valid_aa(2);
is keys %valid_aa, 54;
is($valid_aa{'C'}, 'Cys');
is( $valid_aa{'Cys'}, 'C');


#
# Mutate
#

my $string1 = 'aggt';
$seq = Bio::PrimarySeq->new('-seq'=> 'aggt',
			    '-alphabet'=>'dna',
			    '-id'=>'test3');

# point
Bio::SeqUtils->mutate($seq,
                      Bio::LiveSeq::Mutation->new(-seq => 'c',
                                                  -pos => 3
                                                 )
                     );
is $seq->seq, 'agct';

# insertion and deletion
my @mutations = (
                 Bio::LiveSeq::Mutation->new(-seq => 'tt',
                                             -pos => 2,
                                             -len => 0
                                            ),
                 Bio::LiveSeq::Mutation->new(-pos => 2,
                                             -len => 2
                                            )
);

Bio::SeqUtils->mutate($seq, @mutations);
is $seq->seq, 'agct';

# insertion to the end of the sequence
Bio::SeqUtils->mutate($seq,
                      Bio::LiveSeq::Mutation->new(-seq => 'aa',
                                                  -pos => 5,
                                                  -len => 0
                                                 )
                     );
is $seq->seq, 'agctaa';



#
# testing Bio::SeqUtils->cat
#

# PrimarySeqs

my $primseq1 = Bio::PrimarySeq->new(-id => 1, -seq => 'acgt', -description => 'master');
my $primseq2 = Bio::PrimarySeq->new(-id => 2, -seq => 'tgca');

Bio::SeqUtils->cat($primseq1, $primseq2);
is $primseq1->seq, 'acgttgca';
is $primseq1->description, 'master';

#should work for Bio::LocatableSeq
#should work for Bio::Seq::MetaI Seqs?


# Bio::SeqI

my $seq1 = Bio::Seq->new(-id => 1, -seq => 'aaaa', -description => 'first');
my $seq2 = Bio::Seq->new(-id => 2, -seq => 'tttt', -description => 'second');
my $seq3 = Bio::Seq->new(-id => 3, -seq => 'cccc', -description => 'third');


#  annotations
my $ac2 = Bio::Annotation::Collection->new();
my $simple1 = Bio::Annotation::SimpleValue->new(
                                                -tagname => 'colour',
                                                -value   => 'blue'
                                               ), ;
my $simple2 = Bio::Annotation::SimpleValue->new(
                                                -tagname => 'colour',
                                                -value   => 'black'
                                               ), ;
$ac2->add_Annotation('simple',$simple1);
$ac2->add_Annotation('simple',$simple2);
$seq2->annotation($ac2);

my $ac3 = Bio::Annotation::Collection->new();
my $simple3 = Bio::Annotation::SimpleValue->new(
                                                -tagname => 'colour',
                                                -value   => 'red'
						 ), ;
$ac3->add_Annotation('simple',$simple3);
$seq3->annotation($ac3);


ok (Bio::SeqUtils->cat($seq1, $seq2, $seq3));
is $seq1->seq, 'aaaattttcccc';
is scalar $seq1->annotation->get_Annotations, 3;


# seq features
my $ft2 = Bio::SeqFeature::Generic->new( -start => 1,
                                      -end => 4,
                                      -strand => 1,
                                      -primary => 'source',
                                      -tag     => {note => 'note2'},
				       );


my $ft3 = Bio::SeqFeature::Generic->new( -start => 3,
                                      -end => 3,
                                      -strand => 1,
                                      -primary => 'hotspot',
                                      -tag     => {note => ['note3a','note3b'], 
                                                   comment => 'c1'},
				       );

$seq2->add_SeqFeature($ft2);
$seq2->add_SeqFeature($ft3);


ok (Bio::SeqUtils->cat($seq1, $seq2));
is $seq1->seq, 'aaaattttcccctttt';
is scalar $seq1->annotation->get_Annotations, 5;
is_deeply([uniq_sort(map{$_->get_all_tags}$seq1->get_SeqFeatures)], [sort qw(note comment)], 'cat - has expected tags');
is_deeply([sort map{$_->get_tagset_values('note')}$seq1->get_SeqFeatures], [sort qw(note2 note3a note3b)], 'cat - has expected tag values');
my @tags;
lives_ok {
  @tags = map{$_->get_tag_values(q(note))}$seq1->get_SeqFeatures ;
} 'cat - note tag transfered (no throw)';
cmp_ok(scalar(@tags),'==',3, 'cat - note tag values transfered (correct count)') ;


my $protseq = Bio::PrimarySeq->new(-id => 2, -seq => 'MVTF'); # protein seq

throws_ok {
    Bio::SeqUtils->cat($seq1, $protseq);
} qr/different alphabets:/, 'different alphabets' ;


#
# evolve()
#

$seq = Bio::PrimarySeq->new('-seq'=> 'aaaaaaaaaa',
                            '-id'=>'test');



$util = Bio::SeqUtils->new(-verbose => 0);
ok my $newseq = $util->evolve($seq, 60, 4);

#  annotations

$seq2 = Bio::Seq->new(-id => 2, -seq => 'ggttaaaa', -description => 'second');
$ac3 = Bio::Annotation::Collection->new();
$simple3 = Bio::Annotation::SimpleValue->new(
                                                -tagname => 'colour',
                                                -value   => 'red'
                                                 ), ;
$ac3->add_Annotation('simple',$simple3);
$seq2->annotation($ac3);
$ft2 = Bio::SeqFeature::Generic->new( -start => 1,
                                      -end => 4,
                                      -strand => 1,
                                      -primary => 'source',
                                      -tag     => {note => 'note2'},
                                       );


$ft3 = Bio::SeqFeature::Generic->new( -start => 5,
                                      -end => 8,
                                      -strand => -1,
                                      -primary => 'hotspot',
                                      -tag     => {note => ['note3a','note3b'], 
                                                   comment => 'c1'},
                                       );
$seq2->add_SeqFeature($ft2);
$seq2->add_SeqFeature($ft3);

my $trunc=Bio::SeqUtils->trunc_with_features($seq2, 2, 7);
is $trunc->seq, 'gttaaa';
my @feat=$trunc->get_SeqFeatures;
is $feat[0]->location->to_FTstring, '<1..3';
is $feat[1]->location->to_FTstring, 'complement(4..>6)';
is_deeply([uniq_sort(map{$_->get_all_tags}$trunc->get_SeqFeatures)], [sort qw(note comment)], 'trunc_with_features - has expected tags');
is_deeply([sort map{$_->get_tagset_values('note')}$trunc->get_SeqFeatures], [sort qw(note2 note3a note3b)], 'trunc_with_features - has expected tag values');

my $revcom=Bio::SeqUtils->revcom_with_features($seq2);
is $revcom->seq, 'ttttaacc';
my ($rf1) = $revcom->get_SeqFeatures('hotspot');
is $rf1->primary_tag, $ft3->primary_tag, 'primary_tag matches original feature...';
is $rf1->location->to_FTstring, '1..4', 'but tagged sf is now revcom';

my ($rf2) = $revcom->get_SeqFeatures('source');
is $rf2->primary_tag, $ft2->primary_tag, 'primary_tag matches original feature...';
is $rf2->location->to_FTstring, 'complement(5..8)', 'but tagged sf is now revcom';

is_deeply([uniq_sort(map{$_->get_all_tags}$revcom->get_SeqFeatures)], [sort qw(note comment)], 'revcom_with_features - has expected tags');
is_deeply([sort map{$_->get_tagset_values('note')}$revcom->get_SeqFeatures], [sort qw(note2 note3a note3b)], 'revcom_with_features - has expected tag values');
# check circularity
isnt($revcom->is_circular, 1, 'still not circular');
$seq3 = Bio::Seq->new(-id => 3, -seq => 'ggttaaaa', -description => 'third', -is_circular => 1);
is(Bio::SeqUtils->revcom_with_features($seq3)->is_circular, 1, 'still circular');

sub uniq_sort {
    my @args = @_;
    my %uniq;
    @args = sort @args;
    @uniq{@args} = (0..$#args);
    return sort {$uniq{$a} <=> $uniq{$b}} keys %uniq;
}
