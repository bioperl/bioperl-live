# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }

    use Test;
    plan tests => 12; 


}

if( $error == 1 ) {
    exit(0);
}

use Bio::Tools::tRNAscanSE;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $parser = new Bio::Tools::tRNAscanSE(-verbose => $verbose,
					 -file => Bio::Root::IO->catfile
					 ('t','data', 
					  'yeast.tRNAscanSE'));

my @genes;
while( my $gene = $parser->next_prediction ) {
    push @genes, $gene;
}

ok(@genes, 287);
ok($genes[2]->seq_id, 'I');
my ($codon) = $genes[2]->get_tag_values('Codon');
ok($codon, 'TTG');
ok($genes[2]->start, 181135);
ok($genes[2]->end, 181248);
ok($genes[2]->strand, 1);
ok(my @exons = $genes[2]->get_SeqFeatures,2);
ok($exons[0]->end,181172); 
ok($exons[0]->start,$genes[2]->start); 
ok($exons[1]->start,181205); 
ok($exons[1]->end,$genes[2]->end); 
ok($exons[0]->seq_id, $genes[2]->seq_id);
