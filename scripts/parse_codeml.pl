#!/usr/bin/perl -w
use strict;

use Bio::Tools::Phylo::PAML;
use Bio::Root::IO;

my $parser = new Bio::Tools::Phylo::PAML(-file => Bio::Root::IO->catfile(qw(t data codeml.mlc)));

my $result = $parser->next_result;

my @otus = $result->get_seqs();

my $MLmatrix = $result->get_MLmatrix();
my $NGmatrix = $result->get_NGmatrix();

# These matrices are length(@otu) x length(@otu) "strict lower
# triangle" 2D-matrices, which means that the diagonal and
# everything above it is undefined.  Each of the defined cells is a
# hashref of estimates for "dN", "dS", "omega" (dN/dS ratio), "t",
# "S" and "N".  If a ML matrix, "lnL" will also be defined.
printf "The ML omega ratio for sequences %s vs %s was: %g\n",
    $otus[0]->id, $otus[1]->id, $MLmatrix->[0]->[1]->{omega};

printf "The NG omega ratio for sequences %s vs %s was: %g\n",
    $otus[0]->id, $otus[1]->id, $NGmatrix->[0]->[1]->{'omega'};
