#!/usr/bin/perl
# Author: Jason Stajich <jason.stajich@duke.edu>
# Convert a PAUP tree block to Phylip format

use strict;

my @data;
while(<>) {     
    last if( /Translate/ );
}
while(<>) { 
    last if (/;/);
    my ($num, $taxon) = (/\s+(\d+)\s([A-Za-z\.\_]+),/);
    $data[$num] = substr($taxon,0,10);
}
while(<>) {    
    next unless (s/^\s*tree (\S+) = \[\S+\] //i);
    my $tree = $_;    
    for( my $i=scalar @data; $i > 0; $i-- ) {
	my $taxon = $data[$i];
	$tree =~ s/$i/$taxon/;
    }    
    print $tree;    
}
