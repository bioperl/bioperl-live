#!/bin/perl
# Getting  Entry, Chain, Residue, and Atom objects given a PDB file

use Bio::Structure::IO;
use strict;

my $file = shift or die "No PDB file\n";
my $structio = Bio::Structure::IO->new(-file => $file);
my $struc = $structio->next_structure;

for my $chain ($struc->get_chains) {
   my $chainid = $chain->id;
   # one-letter chaincode if present, 'default' otherwise
   for my $res ($struc->get_residues($chain)) {
      my $resid = $res->id;
      # format is 3-lettercode - dash - residue number, e.g. PHE-20
      my $atoms = $struc->get_atoms($res);
      # actually a list of atom objects, used here to get a count
      print join "\t", $chainid,$resid,$atoms,"\n";
   }
}
