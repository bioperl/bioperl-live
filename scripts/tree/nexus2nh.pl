#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;

my $in = Bio::TreeIO->new(-format => 'nexus', -fh => \*ARGV);
my $out= Bio::TreeIO->new(-format => 'newick');

while( my $t = $in->next_tree ) {
    $out->write_tree($t);
}
