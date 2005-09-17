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
    $NTESTS = 6;
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
my @rpt;
while (my $rpt = $parser->next_result){
    push @rpt, $rpt;
}
ok ($rpt[0]->feature1->start, 1337);
ok ($rpt[0]->feature1->end, 1407);
ok ($rpt[0]->feature1->strand, 1);
ok ($rpt[0]->feature1->primary_tag, "Simple_repeat");
ok ($rpt[0]->feature1->source_tag, "RepeatMasker");
ok ($rpt[0]->feature2->seq_id, "(TTAGGG)n");






