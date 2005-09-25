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
    $NTESTS = 20;
    plan tests => $NTESTS;
}
use Bio::Tools::Pseudowise;
use Bio::Root::IO;

END {
    for ( $Test::ntest..$NTESTS ) {
	#skip("Cannot run remaining pseudowise tests, skipping.",1);
    }
}

my $inputfilename= Bio::Root::IO->catfile("t","data","pseudowise.out");
my $parser = Bio::Tools::Pseudowise->new(-file => $inputfilename);
my @gene;
while (my $gene= $parser->next_prediction){
    push @gene, $gene;
}
my ($g) = @gene;
my @e = $g->sub_SeqFeature;

ok ($g->primary_tag, 'pseudogene');
ok ($g->source_tag, 'pseudowise');

ok(($g->get_tag_values('Synonymous'))[0],7);
ok(($g->get_tag_values('Nonsynonymous'))[0],18);
ok(($g->get_tag_values('Ka/Ks'))[0],2.57);
ok(($g->get_tag_values('Unlikely'))[0],0);
ok(($g->get_tag_values('Identical'))[0],5);
ok(($g->get_tag_values('Stop'))[0],0);
ok(($g->get_tag_values('Total codons'))[0],30);
ok(($g->get_tag_values('Frameshift'))[0],0);
ok(($g->get_tag_values('Intron'))[0],1);

ok($g->start,163);
ok($g->end,626);
ok($g->strand,1);
ok($e[0]->start, 163);
ok($e[0]->end,213);
ok($e[0]->strand,1);
ok($e[1]->start,585);
ok($e[1]->end,626);
ok($e[1]->strand,1);
