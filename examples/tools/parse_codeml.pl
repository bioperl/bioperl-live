#!/usr/bin/perl

use strict;
use Bio::Tools::Phylo::PAML;
use Bio::Root::IO;

my $parser = new Bio::Tools::Phylo::PAML(-file    => shift,
					 -verbose => shift);

my $result = $parser->next_result;
my @otus = $result->get_seqs();
my $MLmatrix = $result->get_MLmatrix();
my $NGmatrix = $result->get_NGmatrix();

# These matrices are length(@otu) x length(@otu) "strict lower
# triangle" 2D-matrices, which means that the diagonal and
# everything above it is undefined.  Each of the defined cells is a
# hashref of estimates for "dN", "dS", "omega" (dN/dS ratio), "t",
# "S" and "N".  If a ML matrix, "lnL" will also be defined.

@otus = $result->get_seqs();
$MLmatrix = $result->get_MLmatrix();
$NGmatrix = $result->get_NGmatrix();
for( my $i=0;$i<scalar @$MLmatrix;$i++) {
	for( my $j = $i+1; $j < scalar @{$MLmatrix->[$i]}; $j++ ) { 
		printf "The ML omega ratio for sequences %s vs %s was: %g\n",
		  $otus[$i]->id, $otus[$j]->id, $MLmatrix->[$i]->[$j]->{omega};
	}
}

for( my $i=0;$i<scalar @$MLmatrix;$i++) {
	for( my $j = $i+1; $j < scalar @{$MLmatrix->[$i]}; $j++ ) { 
	
		printf "The NG omega ratio for sequences %s vs %s was: %g\n",
		  $otus[$i]->id, $otus[$j]->id, $NGmatrix->[$i]->[$j]->{'omega'};
	}
}
