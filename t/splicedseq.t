# -*-Perl-*-

use strict;
use vars qw($DEBUG $TESTCOUNT);
my $error;

BEGIN {     
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 15;
	use_ok('Bio::Seq');
	use_ok('Bio::SeqIO');
};

my $str = Bio::SeqIO->new(
			'-file'=> Bio::Root::IO->catfile("t","data", "U58726.gb"), 
			'-format' => 'GenBank');
ok $str;
my $seq;

ok ( $seq = $str->next_seq() );

# Here is a cute way to verify the sequence by seeing if the
# the translation matches what is annotated in the file -js
foreach my $ft ( grep { $_->primary_tag eq 'CDS'} 
		 $seq->top_SeqFeatures ) {
    if( $ft->has_tag('translation') ) {
	my ($translation) = $ft->each_tag_value('translation');
	my $t = $ft->spliced_seq(-nosort => 1);
	my $pepseq = $t->translate()->seq();
	chop($pepseq);# chop is to remove stop codon
	is($translation,$pepseq); 
	}
}

my $stream = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
                                       ("t","data","M12730.gb"),
                              -format => 'genbank');
# Jump down to M12730 which lists CDS join(1959..2355,1..92)
while ($seq->accession ne "M12730") {
    $seq = $stream->next_seq;
}
ok(my @features = $seq->get_SeqFeatures(), "get_SeqFeatures()");
my $feat;
foreach my $feat2 ( @features ) {
    next unless ($feat2->primary_tag eq "CDS");
    my @db_xrefs = $feat2->annotation->get_Annotations("db_xref");
    if (grep { $_ eq "GI:150830" } @db_xrefs) {
       $feat = $feat2;
       last;
    }
}
my ($protein_seq) = $feat->annotation->get_Annotations("translation");
like($protein_seq, qr(^MKERYGTVYKGSQRLIDE.*ANEKQENALYLIIILSRTSIT$),
	 "protein sequence");
my ($nucleotide_seq) = $feat->spliced_seq(-nosort => 1)->seq;
like($nucleotide_seq, qr(^ATGAAAGAAAGATATGGA.*TCAAGGACTAGTATAACATAA$),
	 "nucleotide sequence - correct CDS range");
is(length($nucleotide_seq), 489, "nucleotide length");
