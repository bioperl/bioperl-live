#!/usr/bin/perl -w
use strict;
use Bio::DB::Taxonomy;

use strict;
use Getopt::Long;
my $verbose = 0;
my $plain   = 0;
my ($nodesfile,$namesfile);

GetOptions('v|verbose' => \$verbose,
	   'nodes:s'   => \$nodesfile,
	   'names:s'   => \$namesfile,
	   'h|help'    => sub{ exec('perldoc',$0);
				exit(0)
				} );

unless( @ARGV || $nodesfile || $namesfile ) {
    exec('perldoc',$0);
    exit(0);
}
my $db = new Bio::DB::Taxonomy(-source    => 'flatfile',
			       -nodesfile => $nodesfile,
			       -namesfile => $namesfile,
			       -directory => '/tmp/idx');
foreach my $sp ( @ARGV ) {
    my $id = $db->get_taxonid($sp);
    print "id is $id for $sp\n";
}

=head1 NAME

local_taxonomydb_query - query a local TaxonomyDB for species or taxonid

=head1 DESCRIPTION

This script provides an example implementation of access to a local
Taxonomy database implemented with Berkeley DB (DB_File module is needed).

Usage:

 local_taxonomydb_query.PLS: [-v] --nodes nodes.dmp --names names.dmp "Genus1 species1" "Genus2 species2"

Providing the nodes.dmp and names.dmp files from the NCBI Taxonomy
dump (see Bio::DB::Taxonomy::flatfile for more info) is only necessary
on the first time running.  This will create the local indexes and may
take quite a long time.  However once created, these indexes will
allow fast access for species to taxon id OR taxon id to species name
lookups.
