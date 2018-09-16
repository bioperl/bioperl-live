#!/usr/bin/perl

# allenday@ucla.edu
# parses a dbsnp xml file, prints some info for each refsnp and subsnp

use strict;
use Bio::ClusterIO;
use Bio::Root::IO;
use IO::File;

my $file = shift @ARGV;

my $io = Bio::ClusterIO->new	(	-tempfile => 0,
					-format   => 'dbsnp',
					-fh       => IO::File->new("zcat $file |"),
				);

while(my $cluster = $io->next_cluster){
	print $cluster->id,"\t", $cluster->observed, "\n";

	foreach my $subsnp ($cluster->each_subsnp){
		print "\t\t\t", $subsnp->id, "\t", $subsnp->handle, "\t", $subsnp->method, "\n";

	}
}
