## Bioperl Test Harness Script for Modules
##


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
BEGIN { $| = 1; print "1..4\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::SeqFeatureProducer::MZEF;
use Bio::SeqFeatureProducer::GFF;
use Bio::SeqIO;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $seqio = new Bio::SeqIO(-format=>'Fasta', -file=>'t/HUMBETGLOA.fasta');
my $seq = $seqio->next_seq;

# trouble with references, this is just easier for now.
$seqio = new Bio::SeqIO(-format=>'Fasta', -file=>'t/HUMBETGLOA.fasta');
my $seq2 = $seqio->next_seq;

# for debugging
$seqio = new Bio::SeqIO(-format=>'GenBank', -fh=>\*STDERR);

print "not ok 2\n" unless ( defined $seq );

print "ok 2\n";
	
my $seqprod = new Bio::SeqFeatureProducer::GFF(-gff=>'t/HUMBETGLOA.gff');
if( ! $seqprod->add_features($seq) ) {
    print "not ok 3\n";
}
else { 
    print "GFF features are for ", $seq->display_id(), "\n";
    foreach my $feat ( $seq->top_SeqFeatures() ) {
	print "Feature from ", $feat->start(), "to ", 
	$feat->end(), " Primary tag  ", $feat->primary_tag(),
	" From", $feat->source_tag(), "\n";
	if( $feat->strand == 0 ) {
	    print "Feature applicable to either strand\n";
	} else {
	    print "Feature on strand ", $feat->strand(),"\n"; # -1,1
	}
	
	foreach $tag ( $feat->all_tags() ) {
	    print "Feature has tag ",$tag,"with value,", $feat->has_tag($tag), "\n";
	}
    }	
}
$seqio->write_seq($seq);
print "ok 3\n";

$seqprod = new Bio::SeqFeatureProducer::MZEF(-mzef=>'t/HUMBETGLOA.mzef');

if( ! $seqprod->add_features($seq2) ) {
    print "not ok 4\n";
}
else { 
    print "\n\nMZEF features are for ", $seq2->display_id(), "\n";
    foreach my $feat ( $seq2->top_SeqFeatures() ) {
	print "Feature from ", $feat->start(), "to ", 
	$feat->end(), " Primary tag  ", $feat->primary_tag(),
	" From", $feat->source_tag(), "\n";
	if( $feat->strand == 0 ) {
	    print "Feature applicable to either strand\n";
	} else {
	    print "Feature on strand ", $feat->strand(),"\n"; # -1,1
	}

	foreach $tag ( $feat->all_tags() ) {
	    print "Feature has tag ",$tag,"with value,", $feat->has_tag($tag), "\n";
	}
    }	
}
$seqio->write_seq($seq2);


print "ok 4\n";
