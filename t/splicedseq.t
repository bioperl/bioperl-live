# -*-Perl-*-

use strict;
use vars qw($DEBUG $TESTCOUNT);
my $error;

BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $TESTCOUNT = 9;
    plan tests => $TESTCOUNT;
    $error = 0;
};

if( $error ==  1 ) {
    exit(0);
}


use Bio::Seq;
use Bio::SeqIO;


my $str = Bio::SeqIO->new( '-file'=> Bio::Root::IO->catfile("t","data",
							    "U58726.gb"), 
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
	my $t = $ft->spliced_seq();
	my $pepseq = $t->translate()->seq();
	chop($pepseq);# chop is to remove stop codon
	ok($translation,$pepseq); 
    }	
}

eval { require Bio::DB::GenBank };
if( $@ ) {
    print STDERR "Skipping remote location tests\n";
    for( $Test::ntest..$TESTCOUNT ) {
	skip("Not possible to test remote locations without DB access",1);
    }
    exit(0);
} else { 
	
#my $db = Bio::DB::GenBank->new();
#
#foreach my $ft ( $seq->top_SeqFeatures ) {
#	my $t = $ft->spliced_seq();
#	print "Got ",$t->seq,"\n";
#}

}


