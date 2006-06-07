# -*-Perl-*-
## Bioperl Test Harness Script for Modules


use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 13;
    plan tests => $NTESTS;
}
use Bio::Tools::RepeatMasker;
use Bio::SeqIO;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("Cannot complete RepeatMasker tests",1);
    }
}

my $inputfilename= Bio::Root::IO->catfile("t","data","repeatmasker.fa.out");
my $parser = Bio::Tools::RepeatMasker->new(-file => $inputfilename);
my $i = 0;
while (my $rpt = $parser->next_result){
    unless( $i++ ) {
	ok ($rpt->feature1->seq_id, "contig11600");
	ok ($rpt->feature1->start, 1337);
	ok ($rpt->feature1->end, 1407);
	ok ($rpt->feature1->strand, 1);
	ok ($rpt->feature1->primary_tag, "Simple_repeat");
	ok ($rpt->feature1->source_tag, "RepeatMasker");
	ok (scalar $rpt->feature1->get_tag_values('Target'), 3);

	ok ($rpt->feature2->seq_id, "(TTAGGG)n");
	ok ($rpt->feature2->start, 2);
	ok ($rpt->feature2->end, 76);
	ok ($rpt->feature1->primary_tag, "Simple_repeat");
	ok ($rpt->feature1->source_tag, "RepeatMasker");
	ok (scalar $rpt->feature2->get_tag_values('Target'), 3);
    }
}






