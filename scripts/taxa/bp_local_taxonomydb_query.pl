#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Taxonomy;

use strict;
use Getopt::Long;
my $verbose = 0;
my $plain   = 0;
my ($nodesfile,$namesfile);
my $idx_dir = '/tmp/idx';
GetOptions('v|verbose' => \$verbose,
	   'nodes:s'   => \$nodesfile,
	   'names:s'   => \$namesfile,
	   'idx:s'     => \$idx_dir,
	   'h|help'    => sub{ exec('perldoc',$0);
				exit(0)
				} );

unless( @ARGV || $nodesfile || $namesfile ) {
    exec('perldoc',$0);
    exit(0);
}
mkdir($idx_dir) unless -d $idx_dir;

my $db = new Bio::DB::Taxonomy(-source    => 'flatfile',
			       -nodesfile => $nodesfile,
			       -namesfile => $namesfile,
			       -directory => $idx_dir);
foreach my $sp ( @ARGV ) {
    my $node = $db->get_Taxonomy_Node(-name => $sp);
    if( defined $node ) {
	print "id is ", $node->id, "\n"; # 9606
	print "rank is ", $node->rank, "\n"; # species
	print "scientific name is ", $node->scientific_name, "\n"; # Homo sapiens
	print "division is ", $node->division, "\n"; # Primates
    } else {
	warn("no node found for query $sp");
    }
}

=head1 NAME

bp_local_taxonomydb_query - query a local TaxonomyDB for species or taxonid

=head1 DESCRIPTION

This script provides an example implementation of access to a local
Taxonomy database implemented with Berkeley DB (DB_File module is needed).

Usage:

 bp_local_taxonomydb_query.PLS: [-v] --nodes nodes.dmp --names names.dmp "Genus1 species1" "Genus2 species2"

Providing the nodes.dmp and names.dmp files from the NCBI Taxonomy
dump (see Bio::DB::Taxonomy::flatfile for more info) is only necessary
on the first time running.  This will create the local indexes and may
take quite a long time.  However once created, these indexes will
allow fast access for species to taxon id OR taxon id to species name
lookups.
