## $Id$

# test for Bio::Seq::PrimedSeq
# written by Rob Edwards

use strict;
use constant NUMTESTS => 9;

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't','..';
    }
    use Test;

    plan tests => NUMTESTS;
}

use Bio::SeqIO;
use Bio::Seq::PrimedSeq;
ok(1);

my ($seqio, $seq, $left, $right, $primed_seq, $left_test, $annseq, $amplicon, $returnedseq);


$seqio=Bio::SeqIO->new(-file=>'t/data/primedseq.fa');
$seq=$seqio->next_seq;
$left=Bio::SeqFeature::Primer->new(-seq=>'CTTTTCATTCTGACTGCAACG');
$right=Bio::SeqFeature::Primer->new(-seq=>'GGTGGTGCTAATGCGT');


ok $primed_seq = Bio::Seq::PrimedSeq->new(-seq=>$seq, -left_primer=>$left, -right_primer=>$right);
ok $left_test = $primed_seq->get_primer('left');
ok $left_test eq $left;
ok $annseq = $primed_seq->annotated_sequence; # should I check that this is what I think it is, or just be happy?
ok $amplicon=$primed_seq->amplicon->seq;
ok uc($amplicon) eq uc('cttttcattctgactgcaacgGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAacgcattagcaccacc');
ok $returnedseq=$primed_seq->seq;
ok $returnedseq->seq eq $seq->seq;

