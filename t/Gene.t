## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..11\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::AnnSeq;
use Bio::Seq;
use Bio::SeqFeature::Exon;
use Bio::SeqFeature::Transcript;
use Bio::SeqFeature::Gene;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $exon = Bio::SeqFeature::Exon->new( -start => 10, -end => 30, -strand => 1);

if( $exon->start == 10 && $exon->end == 30 && $exon->strand == 1 ) {
  print "ok 2\n"; 	
} else {
  print "not ok 2\n"; 	
}


my $trans = Bio::SeqFeature::Transcript->new();

# build up a set of exons into an array
foreach $exonpair ( split(/:/,"10,30:50,80:110,130") ) {
    ($start,$end) = split(/,/,$exonpair);
    push(@exon,Bio::SeqFeature::Exon->new(-start => $start,-end => $end, -strand => 1));
}
    
# add them to a new Transcript
$trans = Bio::SeqFeature::Transcript->new();
$trans->add_Exon(@exon);	

print "ok 3\n";    

($exon1,$exon2,$exon3) = $trans->each_Exon();
#print STDERR "$exon1:$exon2:$exon3 ", join(' ',$exon1->start,$exon2->start,$exon3->start);
if( $exon1->start == 10 && $exon2->start == 50 && $exon3->start == 110 && $trans->start == 10 && $trans->end == 130) {
    print "ok 4\n";
} else {
    print "not ok 4\n";
}

$tls = $trans->create_Translation(20,120);
$trans->Translation($tls);

print "ok 5\n";

#alternative transcript with different translation
pop(@exon);
push(@exon,Bio::SeqFeature::Exon->new(-start => '90',-end => '120',-strand => 1));
# add them to a new Transcript
$alt1 = Bio::SeqFeature::Transcript->new();
$alt1->add_Exon(@exon);	
$alt_tls = $alt1->create_Translation(20,100);
$alt1->Translation($alt_tls);


# alternative transcript with a different last exon
# but the same translation as the first Transcript
$alt2 = Bio::SeqFeature::Transcript->new();
pop(@exon);
push(@exon,Bio::SeqFeature::Exon->new(-start => '110',-end => '135',-strand => 1));
$alt2->add_Exon(@exon);	

# alt has the same translation as trans, even though a different
# splicing structure
$alt2->Translation($tls);

print "ok 6\n";


if( $trans->is_protein_coding == 1 ) {
    print "ok 7\n";
} else {
    print "not ok 7\n";
}

$gene = Bio::SeqFeature::Gene->new();
print "ok 8\n";
$gene->add_Transcript($trans,$alt1,$alt2);
print "ok 9\n";

$seq = Bio::Seq->new(-seq => "GGGATTGAGAGTGATCACTCACGCTAACGTCTGCCCTGTTCCTGTATGGTGAGGCCGCACCACAAGCCACCACCGCCGCCGCCTTCTGCGCAACGCCAACCGCCCGCCAAAACGGATCCTTCCCTGCGCCTGCGCAACCAATCTTGGGACCGGACCTTTTTTCTCCGCCCACTACGCATGCGCAAAGCTAGGACAAACTCCCGCCAACACGCAGGCGCCGTAGGTTCACTGCCTACTCCTGCCCGCCATTTCACGTGTTCTCAGAGGCAGGTGGAACTTCTTAATGCGCCTGCGCAAAACTCGCCATTTTACTACACGTGCGGTCAACAAGAGTTCATTGCAAAAAAATTGTTACCTCCTAGCTGCTTGTCTAATACATAGTGTTAATCATGCTTTGCCAAGCGACTTGACTGTAATATTTGCGCGTGGAAGATTAAAAAGATGTTAAACACCCAAGGTAGATTCAAATGTGAATGATTG", -id => "TestID" );

$aseq = Bio::AnnSeq->new();
$aseq->seq($seq);

$aseq->add_SeqFeature($gene);

print "ok 10\n";

foreach $trans ( $gene->each_Transcript() ) {
    $seq = $trans->seq();
    #print STDERR $seq->out_fasta();
}

foreach $tls ( $gene->each_Translation() ) {
    $seq = $tls->seq();
    #print STDERR $seq->out_fasta();
}

print "ok 11\n";
#print STDERR "\n";
foreach $sf ( $aseq->all_SeqFeatures() ) {
    #print STDERR $sf->gff_string() , "\n";
    $sf->gff_string();
}

