
use Test;

BEGIN { plan tests => 4 };

use Bio::Tools::SiRNA;
use Bio::Seq;
use Bio::SeqIO;

# modules compile
print "ok 1\n";

my $input = Bio::SeqIO->new( -file 	=> 't/data/NM_002254.gb',
			     -format 	=> 'Genbank' );

my $seq = $input->next_seq;

my $sirna = Bio::Tools::SiRNA->new( -target 	=> $seq,
				    );

# first test - cds only
 my @pairs = $sirna->design;

 print "CDS only: got ",scalar(@pairs),"\n";

 if( scalar(@pairs) == 65 ) {
     print "ok 2\n";
 } else {
     print "not ok 2\n";
 }



# next test - include 3prime utr
my @feats = $seq->remove_SeqFeatures;

foreach my $feat (@feats) {
    $seq->add_SeqFeature($feat) unless 
	($feat->primary_tag eq 'Target' or $feat->isa('Bio::SeqFeature::SiRNA::Pair'));
}
$sirna->include_3pr(1);

@pairs2 = $sirna->design;
 print "With 3p UTR: got ",scalar(@pairs2),"\n";

if( scalar(@pairs2) == 140 ) {
    print "ok 3\n";
} else {
    print "not ok 3\n";
}

#third test - naked sequence
my $newseq = Bio::Seq->new( -seq => $seq->seq);

$sirna->target($newseq);

@pairs3 = $sirna->design;
 print "Bare sequence: got ",scalar(@pairs3),"\n";

if( scalar(@pairs3) == 142 ) {
    print "ok 4\n";
} else {
    print "not ok 4\n";
}
