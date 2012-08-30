# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14);
	
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
