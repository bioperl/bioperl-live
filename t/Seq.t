# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

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

    plan tests => 38;
}

use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqFeature::Generic;
use Bio::Annotation;
use Bio::Species;
ok(1);

my $seq = Bio::Seq->new(-seq=>'ACTGTGGCGTCAACT',
                        -desc=>'Sample Bio::Seq object',
			-alphabet => 'dna' );
ok $seq;

my $trunc = $seq->trunc(1,4);
ok $trunc->length,  4, 'truncated sequence was not of length 4';

ok $trunc->seq, 'ACTG', 'truncated sequence was not ACTG instead was '. $trunc->seq();

# test ability to get str function
ok  $seq->seq(),  'ACTGTGGCGTCAACT' ;

ok $seq = Bio::Seq->new(-seq=>'actgtggcgtcaact',
		     -desc=>'Sample Bio::Seq object',
		     -display_id => 'something',
		     -accession_number => 'accnum',
		     -alphabet => 'dna' );

ok uc $seq->alphabet, 'DNA' , 'alphabet was ' .$seq->alphabet();

# basic methods

ok $seq->id(), 'something',  "saw ".$seq->id;
ok $seq->accession_number, 'accnum', "saw ". $seq->accession_number ;
ok $seq->subseq(5, 9),  'tggcg', "subseq(5,9) was ". $seq->subseq(5,9);

my $newfeat = Bio::SeqFeature::Generic->new( -start => 10,
					     -end => 12,
					     -primary => 'silly',
					     -source => 'stuff');


$seq->add_SeqFeature($newfeat);
ok $seq->feature_count, 1;

my $species = new Bio::Species
    (-verbose => 1, 
     -classification => [ qw( sapiens Homo Hominidae
			      Catarrhini Primates Eutheria
			      Mammalia Vertebrata Chordata
			      Metazoa Eukaryota )]);
$seq->species($species);
ok $seq->species->binomial, 'Homo sapiens';
$seq->annotation(new Bio::Annotation('-description' => 'desc-here'));
ok $seq->annotation->description, 'desc-here', 
		 'annotation was ' . $seq->annotation();

#
#  translation tests
#

my $trans = $seq->translate();
ok  $trans->seq(), 'TVAST' , 'translated sequence was ' . $trans->seq();

# unambiguous two character codons like 'ACN' and 'GTN' should give out an amino acid
$seq->seq('ACTGTGGCGTCAAC');
$trans = $seq->translate();
ok $trans->seq(), 'TVAST', 'translated sequence was ' . $trans->seq();

$seq->seq('ACTGTGGCGTCAACA');
$trans = $seq->translate();
ok $trans->seq(), 'TVAST', 'translated sequence was ' . $trans->seq();

$seq->seq('ACTGTGGCGTCAACAG');
$trans = $seq->translate();
ok $trans->seq(), 'TVAST', 'translated sequence was ' . $trans->seq();

$seq->seq('ACTGTGGCGTCAACAGT');
$trans = $seq->translate();
ok $trans->seq(), 'TVASTV', 'translated sequence was ' . $trans->seq();

$seq->seq('ACTGTGGCGTCAACAGTA');
$trans = $seq->translate();
ok $trans->seq(), 'TVASTV', 'translated sequence was ' . $trans->seq();

$seq->seq('AC');
ok $seq->translate->seq , 'T', 'translated sequence was ' . $seq->translate->seq();

#difference between the default and full CDS translation

$seq->seq('atgtggtaa');
$trans = $seq->translate();
ok $trans->seq(), 'MW*' , 'translated sequence was ' . $trans->seq();

$seq->seq('atgtggtaa');
$trans = $seq->translate(undef,undef,undef,undef,1);
ok $trans->seq(), 'MW', 'translated sequence was ' . $trans->seq();

#frame 
my $string;
my @frames = (0, 1, 2);
foreach my $frame (@frames) {
    $string .= $seq->translate(undef, undef, $frame)->seq;
    $string .= $seq->revcom->translate(undef, undef, $frame)->seq;
}
ok $string, 'MW*LPHCGYHVVTT';

#Translating with all codon tables using method defaults
$string = '';
my @codontables = qw( 1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23);
foreach my $ct (@codontables) {
    $string .= $seq->translate(undef, undef, undef, $ct)->seq;
}
ok $string, 'MW*MW*MW*MW*MW*MWQMW*MW*MW*MW*MW*MWYMW*MW*MW*MW*MW*';

# CDS translation set to throw an exception for internal stop codons
$seq->seq('atgtggtaataa');
eval {
    $seq->translate(undef, undef, undef, undef, 'CDS' , 'throw');
};
ok $@ ;

$seq->seq('atgtggtaataa');
ok( $seq->translate('J', '-',)->seq, 'MWJJ');

# tests for RichSeq
my $richseq = Bio::Seq::RichSeq->new( -seq => 'atgtggtaataa',
				      -accession_number => 'AC123',
				      -alphabet => 'rna',
				      -molecule => 'mRNA',		
				      -id => 'id1',
				      -dates => [ '2001/1/1' ],
				      -pid => '887821',
				      -keywords => 'JUNK1 JUNK2',
				      -division => 'Fungi',
				      -secondary_accessions => 'AC1152' );
				 
ok ($richseq);
ok ($richseq->seq, 'atgtggtaataa');
ok ($richseq->display_id, 'id1');
ok (($richseq->get_dates)[0], '2001/1/1');
ok (($richseq->get_secondary_accessions)[0], 'AC1152');
ok ($richseq->accession_number, 'AC123');
ok ($richseq->alphabet, 'rna');
ok ($richseq->molecule, 'mRNA');
ok ($richseq->pid, 887821);
ok ($richseq->division, 'Fungi');
ok ($richseq->keywords, 'JUNK1 JUNK2');
$richseq->seq_version('2');
ok ($richseq->seq_version, 2);
