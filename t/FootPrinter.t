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
    $NTESTS = 26;
    plan tests => $NTESTS;
}
use Bio::Tools::FootPrinter;
use Bio::SeqIO;

END {
	for ( $Test::ntest..$NTESTS ) {
		skip("FootPrinter parser failed.",1);
	}
}


my $inputfilename= Bio::Root::IO->catfile("t","data","footprinter.out");
my $parser = Bio::Tools::FootPrinter->new(-file => $inputfilename);
my @sub;
my @species = qw(TETRAODON CHICKEN MOUSE HAMSTER HUMAN PIG);
while (my $feat = $parser->next_feature){
    ok($feat->seq_id, shift @species);
    foreach my $sub ($feat->sub_SeqFeature){
      push @sub,$sub;
    }
}

ok $sub[0]->seq_id, 'TETRAODON-motif1';
ok $sub[0]->start,352;
ok $sub[0]->end,362;
ok $sub[0]->seq->seq,'tacaggatgca';
ok $sub[1]->seq_id, 'TETRAODON-motif2';
ok $sub[1]->start,363;
ok $sub[1]->end,373;
ok $sub[1]->seq->seq,'ccatatttgga';

ok $sub[2]->seq_id, 'CHICKEN-motif1';
ok $sub[2]->start,363;
ok $sub[2]->end,373;
ok $sub[2]->seq->seq,'cacaggatgta';

ok $sub[3]->seq_id, 'CHICKEN-motif2';
ok $sub[3]->start,376;
ok $sub[3]->end,386;
ok $sub[3]->seq->seq,'ccatataagga';

ok $sub[4]->seq_id, 'MOUSE-motif1';
ok $sub[4]->start,234;
ok $sub[4]->end,243;
ok $sub[4]->seq->seq,'cacaggatgt';













