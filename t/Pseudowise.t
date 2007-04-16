# -*-Perl-*-
## Bioperl Test Harness Script for Modules

use strict;
BEGIN {
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    use vars qw($NTESTS);
    $NTESTS = 22;
    plan tests => $NTESTS;
	use_ok('Bio::Tools::Pseudowise');
	use_ok('Bio::Root::IO');
}


my $inputfilename= Bio::Root::IO->catfile("t","data","pseudowise.out");
my $parser = Bio::Tools::Pseudowise->new(-file => $inputfilename);
my @gene;
while (my $gene= $parser->next_prediction){
    push @gene, $gene;
}
my ($g) = @gene;
my @e = $g->sub_SeqFeature;

is ($g->primary_tag, 'pseudogene');
is ($g->source_tag, 'pseudowise');

is(($g->get_tag_values('Synonymous'))[0],7);
is(($g->get_tag_values('Nonsynonymous'))[0],18);
is(($g->get_tag_values('Ka/Ks'))[0],2.57);
is(($g->get_tag_values('Unlikely'))[0],0);
is(($g->get_tag_values('Identical'))[0],5);
is(($g->get_tag_values('Stop'))[0],0);
is(($g->get_tag_values('Total codons'))[0],30);
is(($g->get_tag_values('Frameshift'))[0],0);
is(($g->get_tag_values('Intron'))[0],1);

is($g->start,163);
is($g->end,626);
is($g->strand,1);
is($e[0]->start, 163);
is($e[0]->end,213);
is($e[0]->strand,1);
is($e[1]->start,585);
is($e[1]->end,626);
is($e[1]->strand,1);
