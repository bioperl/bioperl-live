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

    plan tests => 57;
}

use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqFeature::Generic;
use Bio::Species;
use Bio::Annotation::SimpleValue;

ok(1);

ok my $seq = Bio::Seq->new(-seq=>'ACTGTGGCGTCAACT',
                        -desc=>'Sample Bio::Seq object',
			-alphabet => 'dna',
                        -is_circular => 1
                       );

ok $seq->is_circular;
ok not $seq->is_circular(0);
ok not $seq->is_circular;

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

# check IdentifiableI and DescribableI interfaces
ok $seq->isa('Bio::IdentifiableI');
ok $seq->isa('Bio::DescribableI');
# make sure all methods are implemented
ok $seq->authority("bioperl.org"), "bioperl.org";
ok $seq->namespace("t"), "t";
ok $seq->version(0), 0;
ok $seq->lsid_string(), "bioperl.org:t:accnum";
ok $seq->namespace_string(), "t:accnum.0";
ok $seq->description(), 'Sample Bio::Seq object';
ok $seq->display_name(), "something";

# check that feature accession works regardless of lazy things going on
ok scalar($seq->top_SeqFeatures()), 0;
ok scalar($seq->flush_SeqFeatures()), 0;

my $newfeat = Bio::SeqFeature::Generic->new( -start => 10,
					     -end => 12,
					     -primary => 'silly',
					     -source => 'stuff');


$seq->add_SeqFeature($newfeat);
ok $seq->feature_count, 1;

my $species = Bio::Species->new
    (-verbose => 1, 
     -classification => [ qw( sapiens Homo Hominidae
			      Catarrhini Primates Eutheria
			      Mammalia Vertebrata Chordata
			      Metazoa Eukaryota )]);
$seq->species($species);
ok $seq->species->binomial, 'Homo sapiens';
$seq->annotation->add_Annotation('description',
		 Bio::Annotation::SimpleValue->new(-value => 'desc-here'));
my ($descr) = $seq->annotation->get_Annotations('description');
ok $descr->value(), 'desc-here';
ok $descr->tagname(), 'description';

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
ok ($@ =~ /EX/) ;

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
				      -keywords => 'JUNK1;JUNK2',
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
ok ($richseq->keywords, 'JUNK1; JUNK2');
$richseq->seq_version('2');
ok ($richseq->seq_version, 2);

# tests for subtle misbehaviors
$seq = Bio::Seq->new(-primary_id => 'blah', -accession_number => 'foo');
ok ($seq->accession_number, $seq->primary_seq->accession_number);
ok ($seq->primary_id, $seq->primary_seq->primary_id);
$seq->accession_number('blurb');
$seq->primary_id('bar');
ok ($seq->accession_number, $seq->primary_seq->accession_number);
ok ($seq->primary_id, $seq->primary_seq->primary_id);

