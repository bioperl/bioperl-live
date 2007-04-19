# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    
    if( $@ ) {
	    use lib 't/lib';
    }

    use Test::More;
    plan tests => 15; 

}

use_ok ('Bio::Tools::tRNAscanSE');
use_ok ('Bio::Root::IO');

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $parser = new Bio::Tools::tRNAscanSE(-verbose => $verbose,
					 -file => Bio::Root::IO->catfile
					 ('t','data', 
					  'yeast.tRNAscanSE'));

isa_ok($parser, 'Bio::Tools::tRNAscanSE') ;

my @genes;
while( my $gene = $parser->next_prediction ) {
    push @genes, $gene;
}

is (scalar(@genes), 287);
is($genes[2]->seq_id, 'I', 'seq_id');
my ($codon) = $genes[2]->get_tag_values('Codon');
is($codon, 'TTG', 'codon');
is($genes[2]->start, 181135, 'start');
is($genes[2]->end, 181248, 'end');
is($genes[2]->strand, 1, 'strand');

my @exons = $genes[2]->get_SeqFeatures ;
is ( scalar(@exons), 2, 'exons' );
is($exons[0]->end,181172, 'end' ); 
is($exons[0]->start,$genes[2]->start, 'start'); 
is($exons[1]->start,181205, 'start'); 
is($exons[1]->end,$genes[2]->end, 'end'); 
is($exons[0]->seq_id, $genes[2]->seq_id, 'seq_id');
