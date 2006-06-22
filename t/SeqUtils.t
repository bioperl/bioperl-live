# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##$Id$

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 34;
}

use Bio::PrimarySeq;
use Bio::SeqUtils;
use Bio::LiveSeq::Mutation;
ok 1;

my ($seq, $util, $ascii, $ascii_aa, $ascii3);

#                     !    !          
$ascii =    'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
$ascii_aa = 'ABCDEFGHIXKLMNXPQRSTUVWXYZ';

$ascii3 = 
    'AlaAsxCysAspGluPheGlyHisIleXaaLysLeuMetAsnXaaProGlnArgSerThrSecValTrpXaaTyrGlx';

$seq = Bio::PrimarySeq->new('-seq'=> $ascii,
			    '-alphabet'=>'protein', 
			       '-id'=>'test');

# one letter amino acid code to three letter code
ok $util = new Bio::SeqUtils;
ok $util->seq3($seq), $ascii3;

#using anonymous hash
ok (Bio::SeqUtils->seq3($seq), $ascii3); 
ok (Bio::SeqUtils->seq3($seq, undef, ','), 
    'Ala,Asx,Cys,Asp,Glu,Phe,Gly,His,Ile,Xaa,Lys,'.
    'Leu,Met,Asn,Xaa,Pro,Gln,Arg,Ser,Thr,Sec,Val,Trp,Xaa,Tyr,Glx');

$seq->seq('asd-KJJK-');
ok (Bio::SeqUtils->seq3($seq, '-', ':'), 
    'Ala:Ser:Asp:Ter:Lys:Xaa:Xaa:Lys:Ter');

# three letter amino acid code to one letter code
ok (Bio::SeqUtils->seq3in($seq, 'AlaPYHCysAspGlu')), 
ok  $seq->seq, 'AXCDE';
ok (Bio::SeqUtils->seq3in($seq, $ascii3)->seq, $ascii_aa);
#ok ();

#
# Tests for multiframe translations
#

$seq = Bio::PrimarySeq->new('-seq'=> 'agctgctgatcggattgtgatggctggatggcttgggatgctgg',
			    '-alphabet'=>'dna', 
			    '-id'=>'test2');

my @a = $util->translate_3frames($seq);
ok scalar @a, 3;
#foreach $a (@a) {
#    print 'ID: ', $a->id, ' ', $a->seq, "\n";
#}

@a = $util->translate_6frames($seq);
ok scalar @a, 6;
#foreach $a (@a) {
#    print 'ID: ', $a->id, ' ', $a->seq, "\n";
#}

#
# test for valid AA return
#

my @valid_aa = sort Bio::SeqUtils->valid_aa;
ok(@valid_aa, 25);
ok ($valid_aa[1], 'A');

@valid_aa = sort Bio::SeqUtils->valid_aa(1);
ok(@valid_aa, 25);
ok ($valid_aa[1], 'Arg');

my %valid_aa = Bio::SeqUtils->valid_aa(2);
ok keys %valid_aa, 50;
ok($valid_aa{'C'}, 'Cys');
ok( $valid_aa{'Cys'}, 'C');


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
ok $seq->seq, 'agct';

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
ok $seq->seq, 'agct';

# insertion to the end of the sequence
Bio::SeqUtils->mutate($seq,
                      Bio::LiveSeq::Mutation->new(-seq => 'aa',
                                                  -pos => 5,
                                                  -len => 0
                                                 )
                     );
ok $seq->seq, 'agctaa';



#
# testing Bio::SeqUtils->cat
#

use Bio::Annotation::SimpleValue;
use Bio::Seq::RichSeq;;


# PrimarySeqs

my $primseq1 = new Bio::PrimarySeq(-id => 1, -seq => 'acgt', -description => 'master');
my $primseq2 = new Bio::PrimarySeq(-id => 2, -seq => 'tgca');

Bio::SeqUtils->cat($primseq1, $primseq2);
ok $primseq1->seq, 'acgttgca';
ok $primseq1->description, 'master';

#should work for Bio::LocatableSeq
#should work for Bio::Seq::MetaI Seqs?


# Bio::SeqI

my $seq1 = new Bio::Seq(-id => 1, -seq => 'aaaa', -description => 'first');
my $seq2 = new Bio::Seq(-id => 2, -seq => 'tttt', -description => 'second');
my $seq3 = new Bio::Seq(-id => 3, -seq => 'cccc', -description => 'third');


#  annotations
my $ac2 = new Bio::Annotation::Collection;
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

my $ac3 = new Bio::Annotation::Collection;
my $simple3 = Bio::Annotation::SimpleValue->new(
                                                -tagname => 'colour',
                                                -value   => 'red'
						 ), ;
$ac3->add_Annotation('simple',$simple3);
$seq3->annotation($ac3);


ok (Bio::SeqUtils->cat($seq1, $seq2, $seq3));
ok $seq1->seq, 'aaaattttcccc';
ok scalar $seq1->get_Annotations, 3;


# seq features
use Bio::SeqFeature::Generic;

my $ft2 = new Bio::SeqFeature::Generic ( -start => 1,
                                      -end => 4,
                                      -strand => 1,
                                      -primary => 'source',
				       );


my $ft3 = new Bio::SeqFeature::Generic ( -start => 3,
                                      -end => 3,
                                      -strand => 1,
                                      -primary => 'hotspot',
				       );

$seq2->add_SeqFeature($ft2);
$seq2->add_SeqFeature($ft3);


ok (Bio::SeqUtils->cat($seq1, $seq2));
ok $seq1->seq, 'aaaattttcccctttt';
ok scalar $seq1->get_Annotations, 5;


my $protseq = new Bio::PrimarySeq(-id => 2, -seq => 'MVTF'); # protein seq

eval {
    Bio::SeqUtils->cat($seq1, $protseq);
};
ok 1 if $@; # did throw

#use Data::Dumper; print Dumper $seq1;






#
# evolve()
#

$seq = Bio::PrimarySeq->new('-seq'=> 'aaaaaaaaaa',
                            '-id'=>'test');



$util = new Bio::SeqUtils(-verbose => 0);
ok my $newseq = $util->evolve($seq, 60, 4);

#  annotations

$seq2 = new Bio::Seq(-id => 2, -seq => 'ttttaaaa', -description => 'second');
$ac3 = new Bio::Annotation::Collection;
$simple3 = Bio::Annotation::SimpleValue->new(
                                                -tagname => 'colour',
                                                -value   => 'red'
                                                 ), ;
$ac3->add_Annotation('simple',$simple3);
$seq2->annotation($ac3);
$ft2 = new Bio::SeqFeature::Generic ( -start => 1,
                                      -end => 4,
                                      -strand => 1,
                                      -primary => 'source',
                                       );


$ft3 = new Bio::SeqFeature::Generic ( -start => 5,
                                      -end => 8,
                                      -strand => -1,
                                      -primary => 'hotspot',
                                       );
$seq2->add_SeqFeature($ft2);
$seq2->add_SeqFeature($ft3);

my $trunc=Bio::SeqUtils->trunc_with_features($seq2, 2, 7);
ok $trunc->seq, 'tttaaa';
my @feat=$trunc->get_SeqFeatures;
ok $feat[0]->location->to_FTstring, '<1..3';
ok $feat[1]->location->to_FTstring, 'complement(4..>6)';
