#!/usr/bin/perl
use strict;

# This program will read in an NMR derived PDB file containing
# multiple conformers, and will allow the user to extract either
# one or all of the models to separate files.
#
# Although the program will run interactively, the command line
# usage is "nmrsplit [input.file] [number of model to extract]"
#
#
#
# Written 13/12/00 by Simon Andrews (simon.andrews@bbsrc.ac.uk)

# Submitted to bioperl script project 2001/08/06

# Description:  Script which reads  an NMR-derived  multiple structure
# PDB file, and will either extract a single structure from it, or
# will  extract all of the structures into single  files.  This  is
# useful when you want to  work with a single representative structure
# from an NMR  ensemble - especially in conjunction  with the OLDERADO
# database (http://neon.chem.le.ac.uk/olderado/) which  finds  the
# most representative structure from an ensemble.


my $Input = $ARGV[0];  	# File to be read
my $Pullout = $ARGV[1];	# Specific model to extract
my @Pdbfile; 		# Array for whole PDB file
my $Header = ""; 	# String to hold the PDB header
my $Model = ""; 	# String to hold individual models
my $Output; 		# Prefix for output files
my $Modno = 1; 		# Number of the model being processed

while () {
    if ($Input)	{

	if (-r $Input) {
	    last;
	}else{
	    print "\"$Input\" does not exist, or could not be read\n";
	}
    }

    print "\nEnter name of multiple PDB file: ";
    $Input = <STDIN>;
    chomp $Input;
    $Input =~ s/^\s*//;
    $Input =~ s/\s*$//;

    next;
}

while () {
    if ($Pullout) {

	if ($Pullout =~ /^\d+$/){
	    if ($Pullout == int $Pullout) {
		last;
	    }else {
		print "\"$Pullout\" should be an integer\n";
	    }
	}else  {
	    print "\"$Pullout\" should be a number\n";
	}
    }

    print "\nEnter number of specific model to extract (Return for none): ";
    $Pullout = <STDIN>;
    chomp $Pullout;
    $Pullout =~ s/^\s*//;
    $Pullout =~ s/\s*$//;

    last unless ($Pullout);
    next;
}


($Output = $Input) =~ s/\.\w*$//; # Take off everything after the last . to use as prefix


open (PDB,$Input) || die "Can't open $Input because $!";



########## Read the header information ####################


while (<PDB>) {
    if (/^MODEL\b/){last;}
    $Header = $Header . $_;
}


######### Read the separate models #######################


while () {

    model();
    if ($Model)	{      # Check if we're past the last model
	if ($Pullout) {		 # Check if we're writing one or all
	    last if ($Modno > $Pullout);# No point continuing if we've got the one we want
	    readout();
	}else	{
	    writeout();
	}
	$Model = "";
	++$Modno;
    }else  {
	last;
    }
}
--$Modno;     # Correct last increment which didn't find a model

if (($Pullout) & ($Modno < $Pullout)) {
    print "\nCannot find model $Pullout : Only $Modno models in this file\n";
}

#################### subroutines start here ##########################


sub model {

    while (<PDB>) {
	if (/^(MODEL\b|END\b|MASTER\b)/){next;} 
	# Stops you getting MODEL... at the top of the output
	# and makes sure there isn't a file containing just END or MASTER

	if (/^ENDMDL\b/){last;}    # Check for the end of the model
	$Model = $Model . $_; 	   # Append the line to $Model
    }
}

sub writeout  {	 # Used when all files are being written out

    if (-e "$Output\_$Modno.pdb"){  # Check whether we're overwriting anything

	print "\n$Output\_$Modno.pdb already exists. Overwrite (y/n)? ";
	my $Question = <STDIN>;
	unless ($Question =~ /^y/i) {
	    print "\nSkipping $Output\_$Modno.pdb";
	    return;
	}
    }

    open (OUT,">$Output\_$Modno.pdb") || die "Can't open $Output\_$Modno.pdb because $!";
    print "\nWriting $Output\_$Modno.pdb ...";
    print OUT $Header;
    print OUT $Model;
    print OUT "END\n";	# Adds and END statement to the PDB file

    close OUT || die "Couldn't close $Output\_$Modno.pdb because $!";
}


sub readout {

    if ($Modno == $Pullout) {

	if (-e "$Output\_$Modno.pdb") {	# Check whether we're overwriting anything

	    print "\n$Output\_$Modno.pdb already exists. Overwrite (y/n)? ";
	    my $Question = <STDIN>;
	    unless ($Question =~ /^y/i)	{
		print "\nModel not extracted\n";
		$Model = "";
		return;
	    }
	}

	open (OUT,">$Output\_$Modno.pdb") || die "Can't open $Output\_$Modno.pdb because $!";
	print "\nWriting $Output\_$Modno.pdb ...\n";
	print OUT $Header;
	print OUT $Model;
	print OUT "END\n"; # Adds and END statement to the PDB file

	close OUT || die "Couldn't close $Output\_$Modno.pdb because $!";

	$Model = "";  # Stops the reading after this model
    }else  {
	print "\nReading Model $Modno ...";
    }
}
