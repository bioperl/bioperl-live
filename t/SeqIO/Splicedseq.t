# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 15);
	
	use_ok('Bio::SeqIO');
}

ok my $str = Bio::SeqIO->new(
			'-file'=> test_input_file('U58726.gb'), 
			'-format' => 'GenBank');

my $seq;

ok ( $seq = $str->next_seq() );

# Here is a cute way to verify the sequence by seeing if the
# the translation matches what is annotated in the file -js
foreach my $ft ( grep { $_->primary_tag eq 'CDS'} 
		 $seq->top_SeqFeatures ) {
    if( $ft->has_tag('translation') ) {
	my ($translation) = $ft->get_tag_values('translation');
	my $t = $ft->spliced_seq(-nosort => 1);
	my $pepseq = $t->translate()->seq();
	chop($pepseq);# chop is to remove stop codon
	is($translation,$pepseq); 
	}
}

my $stream = Bio::SeqIO->new(-file => test_input_file('M12730.gb'),
                              -format => 'genbank');
# Jump down to M12730 which lists CDS join(1959..2355,1..92)
while ($seq->accession ne "M12730") {
    $seq = $stream->next_seq;
}
ok(my @features = $seq->get_SeqFeatures(), "get_SeqFeatures()");
my $feat;
foreach my $feat2 ( @features ) {
    next unless ($feat2->primary_tag eq "CDS");
    my @db_xrefs = $feat2->get_tag_values("db_xref");
    if (grep { $_ eq "GI:150830" } @db_xrefs) {
       $feat = $feat2;
       last;
    }
}
my ($protein_seq) = $feat->get_tag_values("translation");
like($protein_seq, qr(^MKERYGTVYKGSQRLIDE.*ANEKQENALYLIIILSRTSIT$),
	 "protein sequence");
my ($nucleotide_seq) = $feat->spliced_seq(-nosort => 1)->seq;
like($nucleotide_seq, qr(^ATGAAAGAAAGATATGGA.*TCAAGGACTAGTATAACATAA$),
	 "nucleotide sequence - correct CDS range");
is(length($nucleotide_seq), 489, "nucleotide length");

# Test remote sequence handling
my $db_in;
eval{
    use Bio::DB::GenBank;
    use Test::Warn;
    $db_in = Bio::DB::GenBank -> new;
    my $seq_obj = $db_in->get_Seq_by_id('AF032048.1');
};
if($@){
    warn("Warning: Problem with accessing GenBank entry AF032048.1 to test spliced_seq on remote DBs. skipping this test\n");
} else {
    my $str = Bio::SeqIO->new(
			'-file'=> test_input_file('AF032047.gbk'), 
			'-format' => 'GenBank');
    my @feats = $str-> next_seq -> get_SeqFeatures;
    # feat[1] has some exons from external sequence AF032048.1
    # Compare output seq length with and without specifying db (latter results in N padding)
    
    my $len_w_db = length($feats[1]->spliced_seq(-db => $db_in)->seq);
    my $len_nodb = length($feats[1]->spliced_seq(             )->seq);     # produces 2 warnings
    ok($len_nodb == $len_w_db, "right number of Ns added if source sequence not found");
}
