#!/usr/bin/perl -w

# set this if BioPerl is not installed in a standard location
#use lib "/home/kris/projects/bioperl/HEAD/bioperl-live";

# set this to the location of the PDB file. An example can be found
# in the BioPerl distribution in t/data
my $pdb_file = "../../t/data/pdb1bpt.ent";

#
# nothing needs to be changed below here
#

# This simple script reads in a PDB file and loops over all residues
# in it. The atom name and the x,y,z coordinates are then printed.
# (similar methods to loop over all chains exist as well)
#
# Author: Kris Boulez (kris.boulez@algonomics.com)
# $Id$

use strict;

use Bio::Structure::IO;

# create a structIO handle
my $io = new Bio::Structure::IO->new(-file => $pdb_file,
                                     -format => 'pdb');

# read in the first structure
my $struc = $io->next_structure;

my @residues = $struc->residue;
print "residue count is ", scalar @residues, "\n";

# loop over all residues
for my $residue (@residues) {
        print "residue: ", $residue->id,"\n";
        # get all atoms for this residue
        for my $atom ( $struc->get_atoms($residue) ) {
                print "\tatom id is ", $atom->id, " x is ", $atom->x,
                   " y is ", $atom->y, " z is ", $atom->z, "\n";
        }
} 

