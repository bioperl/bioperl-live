# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 73);

    use_ok('Bio::Seq');
    use_ok('Bio::Seq::RichSeq');
    use_ok('Bio::SeqFeature::Generic');
    use_ok('Bio::Species');
    use_ok('Bio::Annotation::SimpleValue');
}

ok my $seq = Bio::Seq->new(-seq=>'ACTGTGGCGTCAACT',
                        -desc=>'Sample Bio::Seq object',
                        -alphabet => 'dna',
                        -is_circular => 1
                       );
isa_ok($seq,"Bio::AnnotatableI");

ok $seq->is_circular;
ok not $seq->is_circular(0);
ok not $seq->is_circular;

my $trunc = $seq->trunc(1,4);
is $trunc->length, 4, 'truncated sequence length';

is $trunc->seq, 'ACTG', 'truncated sequence string';

# test ability to get str function
is $seq->seq(),  'ACTGTGGCGTCAACT' ;

ok $seq = Bio::Seq->new(-seq=>'actgtggcgtcaact',
                        -desc=>'Sample Bio::Seq object',
                        -display_id => 'something',
                        -accession_number => 'accnum',
                        -alphabet => 'dna' );

is uc $seq->alphabet, 'DNA' , 'alphabet';

# basic methods

is $seq->id(), 'something',  "id";
is $seq->accession_number, 'accnum', "accession number";
is $seq->subseq(5, 9),  'tggcg', "subseq";

# check IdentifiableI and DescribableI interfaces
isa_ok $seq, 'Bio::IdentifiableI';
isa_ok $seq, 'Bio::DescribableI';
# make sure all methods are implemented
is $seq->authority("bioperl.org"), "bioperl.org";
is $seq->namespace("t"), "t";
is $seq->version(0), 0;
is $seq->lsid_string(), "bioperl.org:t:accnum";
is $seq->namespace_string(), "t:accnum.0";
is $seq->description(), 'Sample Bio::Seq object';
is $seq->display_name(), "something";

# check that feature accession works regardless of lazy things going on
is scalar($seq->top_SeqFeatures()), 0;
is scalar($seq->flush_SeqFeatures()), 0;

my $newfeat = Bio::SeqFeature::Generic->new( -start => 10,
                                             -end => 12,
                                             -primary => 'silly',
                                             -source => 'stuff');


$seq->add_SeqFeature($newfeat);
is $seq->feature_count, 1;

my $species = Bio::Species->new
    (-verbose => 1,
     -classification => [ qw( sapiens Homo Hominidae
                              Catarrhini Primates Eutheria
                              Mammalia Vertebrata Chordata
                              Metazoa Eukaryota )]);
$seq->species($species);
is $seq->species->binomial, 'Homo sapiens';
$seq->annotation->add_Annotation('description',
                 Bio::Annotation::SimpleValue->new(-value => 'desc-here'));
my ($descr) = $seq->annotation->get_Annotations('description');
is $descr->value(), 'desc-here';
is $descr->tagname(), 'description';

#
#  translation tests
#

my $trans = $seq->translate();
is  $trans->seq(), 'TVAST' , 'translated sequence';

# unambiguous two character codons like 'ACN' and 'GTN' should give out an amino
# acid ...with the addendum that there should be no assumption by the method
# to complete the codon unless specified, using the -complete_codons flag.

$seq->seq('ACTGTGGCGTCAACN');
$trans = $seq->translate();
is $trans->seq(), 'TVAST', 'translated sequence with explicit unambiguous codons';

$seq->seq('ACTGTGGCGTCAAC');
$trans = $seq->translate();
is $trans->seq(), 'TVAS', 'translated sequence with unknown unambiguous codons';

$seq->seq('ACTGTGGCGTCAAC');
$trans = $seq->translate(-complete_codons => 1);
is $trans->seq(), 'TVAST', 'translated sequence with unknown unambiguous codons, completed';

$seq->seq('ACTGTGGCGTCAACA');
$trans = $seq->translate();
is $trans->seq(), 'TVAST', 'translated sequence with unambiguous codons';

$seq->seq('ACTGTGGCGTCAACAG');
$trans = $seq->translate();
is $trans->seq(), 'TVAST', 'translated sequence with unambiguous codons';

$seq->seq('ACTGTGGCGTCAACAGT');
$trans = $seq->translate(-complete_codons => 1);
is $trans->seq(), 'TVASTV', 'translated sequence with unknown unambiguous codons, completed';

$seq->seq('ACTGTGGCGTCAACAGTA');
$trans = $seq->translate();
is $trans->seq(), 'TVASTV', 'translated sequence with unambiguous codons';

$seq->seq('AC');
is $seq->translate(-complete_codons => 1)->seq , 'T', 'translated sequence with unknown unambiguous codons, completed';

#difference between the default and full CDS translation

$seq->seq('atgtggtaa');
$trans = $seq->translate();
is $trans->seq(), 'MW*' , 'translated sequence with stop';

$seq->seq('atgtggtaa');
$trans = $seq->translate(undef,undef,undef,undef,1);
is $trans->seq(), 'MW', 'translated sequence';

#frame
my $string;
my @frames = (0, 1, 2);
foreach my $frame (@frames) {
    $string .= $seq->translate(undef, undef, $frame)->seq;
    $string .= $seq->revcom->translate(undef, undef, $frame)->seq;
}
is $string, 'MW*LPHCGYHVVTT';

#Translating with all codon tables using method defaults
$string = '';
my @codontables = qw( 1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23);
foreach my $ct (@codontables) {
    $string .= $seq->translate(undef, undef, undef, $ct)->seq;
}
is $string, 'MW*MW*MW*MW*MW*MWQMW*MW*MW*MW*MW*MWYMW*MW*MW*MW*MW*';

# CDS translation set to throw an exception for internal stop codons
$seq->seq('atgtggtaataa');
eval {
    $seq->translate(undef, undef, undef, undef, 'CDS' , 'throw');
};
like ($@, qr/EX/);

$seq->seq('atgtggtaataa');
is( $seq->translate('J', '-',)->seq, 'MWJJ');

# tests for RichSeq
ok my $richseq = Bio::Seq::RichSeq->new( -seq => 'atgtggtaataa',
                                      -accession_number => 'AC123',
                                      -alphabet => 'rna',
                                      -molecule => 'mRNA',
                                      -id => 'id1',
                                      -dates => [ '2001/1/1' ],
                                      -pid => '887821',
                                      -keywords => 'JUNK1;JUNK2',
                                      -division => 'Fungi',
                                      -secondary_accessions => 'AC1152' );

is ($richseq->seq, 'atgtggtaataa');
is ($richseq->display_id, 'id1');
is (($richseq->get_dates)[0], '2001/1/1');
is (($richseq->get_secondary_accessions)[0], 'AC1152');
is ($richseq->accession_number, 'AC123');
is ($richseq->alphabet, 'rna');
is ($richseq->molecule, 'mRNA');
is ($richseq->pid, 887821);
is ($richseq->division, 'Fungi');
is ($richseq->keywords, 'JUNK1; JUNK2');
$richseq->seq_version('2');
is ($richseq->seq_version, 2);

# tests for subtle misbehaviors
$seq = Bio::Seq->new(-primary_id => 'blah', -accession_number => 'foo');
is ($seq->accession_number, $seq->primary_seq->accession_number);
is ($seq->primary_id, $seq->primary_seq->primary_id);
$seq->accession_number('blurb');
$seq->primary_id('bar');
is ($seq->accession_number, $seq->primary_seq->accession_number);
is ($seq->primary_id, $seq->primary_seq->primary_id);


# Bug #2864:

$seq = Bio::Seq->new(-display_id => 0, -seq => 'GATC');

is $seq->display_id, 0, "Bug #2864";

# transcribe/rev_transcribe

$seq = Bio::Seq->new( -id => 'seq1', -alphabet=>'dna',
                      -seq=> 'attTcgcatgT' );
ok my $xseq = $seq->transcribe;
is $xseq->alphabet, 'rna';
ok !($xseq->seq =~ /[tT]/);
is $xseq->seq, 'auuUcgcaugU';
ok !$xseq->transcribe;
ok $seq = $xseq->rev_transcribe;
is $seq->seq, 'attTcgcatgT';
is $seq->alphabet, 'dna';
