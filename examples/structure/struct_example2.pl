#!/usr/local/bin/perl -w

# set this if your BioPerl is not installed in a standard location
# use lib "/home/kris/projects/bioperl/HEAD/bioperl-live";

# the location of a PDB structure file. The file for BPTI
# can be found in the t/data directory of the Bioperl distribution
my $bpti_file = "../../t/data/pdb1bpt.ent";

# nothing needs to be changed below here



#
# We read in the PDB structure for BPTI (1bpt). The structure
# contains 6 cysteines. We check all possible combinations of 
# disulfide bridges by looking at the distance between the SG
# atoms.
# After that we look at the annotation in PDB file
# 
# Author: Kris Boulez (kris.boulez@algonomics.com)
# $Id$

use Bio::Structure::IO;

# create structIO handle
my $structio = Bio::Structure::IO->new(-file => $bpti_file, 
					  -format => 'pdb');

# read the structure
my $struc = $structio->next_structure;

# loop over whole structure and store CYS SG atoms
my @sgatoms;
for my $res ($struc->residue) {
	# we only look at CYS residues
	next unless($res->id =~ /^CYS/);
	# and we only take SG atoms
	for my $atom ($struc->get_atoms($res)) {
		next unless( $atom->id eq "SG");
		push @sgatoms, $atom;
	}
}

print "The measured SG-SG distances are\n--------------------------------\n";

# we now loop over all SG atoms and calculate the distance between them
for my $k (0 .. $#sgatoms) {
	# to enhance readibility of the code below
	my $atom1 = $sgatoms[$k];
	for my $i ($k+1 .. $#sgatoms) {
		my $atom2 = $sgatoms[$i];
		my $dist = calculate_distance($atom1,$atom2);
		printf("%-6s %s - %-6s %s   %-.2f\n",
			$struc->parent($atom1)->id, $atom1->id, 
			$struc->parent($atom2)->id,$atom2->id, $dist);
	}
	print "\n";
}

# now have a look what the annotation whas in the PDB file
print "\nThe annotation in the PDB file say\n";
print "----------------------------------\n";
my ($ann) = $struc->annotation->get_Annotations("ssbond");
if ( defined $ann) {
	# there were SSBOND record in the PDB file
	#
	# the Annotation object has an as_text method to get the text
	my $txt = $ann->as_text;
	# this text starts with "Value: "
	$txt =~ s/^Value: //;
	# and it contains lines of 65 chars long each
	for (my $t = 0; $t <= length($txt); $t += 65) {
		my $line = substr ($txt,$t, 65);
		print "$line\n";
	}
} else {
	print "No annotation found about SSBOND\n";
}


sub calculate_distance {
	# calculate distance between two Atom objects
	my ($atom1, $atom2) = @_;
	my $dist;
	$dist = sqrt( ($atom1->x - $atom2->x)**2 +
		      ($atom1->y - $atom2->y)**2 +
		      ($atom1->z - $atom2->z)**2 );
}
