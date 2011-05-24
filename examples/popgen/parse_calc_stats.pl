#!/usr/bin/perl
# Author: Jason Stajich, jason@bioperl.org
# $Revision: 6576 $

use strict;

use Bio::PopGen::IO;
use Bio::PopGen::Statistics;
use Bio::PopGen::Population;

my $io = new Bio::PopGen::IO(-format => 'prettybase',
			     # the Bio::Root::IO->catfile is only
			     # to make file access platform independent
			     -file   => Bio::Root::IO->catfile
			     (qw( t data popstats.prettybase)));

# This is an example of how to read in data from Bio::PopGen::IO
# We're going to make 2 lists, @outgroup, @ingroup
# @outgroup is a single individual which is named 'out'
# @ingroup is the set of individuals we are testing
my (@ingroup,@outgroup);
while( my $ind = $io->next_individual ) {
    if($ind->unique_id =~ /out/) {
	push @outgroup, $ind;
    } else { 
	push @ingroup, $ind;	
    }
}

# We'll get the names of all the markers (or sites)
# that this individual has genotypes for
my @marker_names = $ingroup[0]->get_marker_names();

# the number of sites is the same as the number of markers
# we assume that all the individuals have the same number of sites
# or that this data is 'aligned' if these were derived from a 
# multiple sequence alignment
my $sitecount = scalar @marker_names;

foreach my $ind ( @ingroup ) {
    # here let's print out the individual name and all their alleles
    # for all the markers
    # like this
    # Name: INDIVIDUALNAME
    #      A1,A2 B1,B2,...
    print "Name: ", $ind->unique_id,"\n";
    print "\t";
    foreach my $marker ( @marker_names ) {
	for my $genotype ( $ind->get_Genotypes($marker) ) {
	    my @alleles = $genotype->get_Alleles();
	    # In this example these are actually single alleles anyways...
	    print join(",", @alleles), " ";
	}
    }
    print "\n";
    
    # There is a more compact way to write that
    print "Name: ", $ind->unique_id,
          "\n\t", join(" ", map { join(",",$_->get_Alleles) } 
		          map { $ind->get_Genotypes($_) } @marker_names),"\n";
    print "--\n";
}

# We can compute some statistics about these individuals
# (underlying assumption is that they are unrelated...)

print "Pi: ",Bio::PopGen::Statistics->pi(\@ingroup), "\n";
print "Theta: ",Bio::PopGen::Statistics->theta(\@ingroup), "\n";

# we can also treat them like a population
my $ingroup_pop = new Bio::PopGen::Population(-individuals => \@ingroup);

print "Pi: ",Bio::PopGen::Statistics->pi($ingroup_pop), "\n";
print "Theta: ",Bio::PopGen::Statistics->theta($ingroup_pop), "\n";





# You can also simulate individuals from a coalescent 
use Bio::PopGen::Simulation::Coalescent;

my $ssize = 5;
my $sim = new Bio::PopGen::Simulation::Coalescent(-sample_size => $ssize);
my $tree = $sim->next_tree;
my $mutcount = 100;
$sim->add_Mutations($tree, $mutcount);

# The leaves are the simulated individuals
my @leaves = $tree->get_leaf_nodes;

# We can use the Stats module either like Bio::PopGen::Statistics->XXX
# or like this:
my $stats = new Bio::PopGen::Statistics;
# $stats->verbose(1);
print "Coalescent pi: ", $stats->pi(\@leaves), "\n";
print "Coalescent theta: ", $stats->theta(\@leaves), "\n";
my $coalescent_pop = new Bio::PopGen::Population(-individuals => \@leaves);

print "Coalescent pi: ", $stats->pi($coalescent_pop), "\n";
print "Coalescent theta: ", $stats->theta($coalescent_pop), "\n";
