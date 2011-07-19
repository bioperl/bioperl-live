#!/usr/bin/perl

=head1 NAME

taxonomy2tree - Building a taxonomic tree based on the full lineages of a set of species names

=head1 DESCRIPTION

taxonomy2tree.PLS: -s Orangutan -s Gorilla -s Chimpanzee -s Human
taxonomy2tree.PLS: -s Orangutan -s Gorilla -s Chimpanzee -s "Homo Sapiens"

Can also provide -d to specific the directory to store index files in, -o to
specific the location of your nodes file, and -a for the names file.
Or the option -e to use the web-based entrez taxonomy database if you don't
have the flatfiles installed.

This script requires that the bioperl-run pkg be also installed.

Providing the nodes.dmp and names.dmp files from the NCBI Taxonomy
dump (see Bio::DB::Taxonomy::flatfile for more info) is only necessary
on the first time running.  This will create the local indexes and may
take quite a long time.  However once created, these indexes will
allow fast access for species to taxon id OR taxon id to species name
lookups.

=head1 AUTHOR - Gabriel Valiente, reimplemented by Sendu Bala

Email valiente@lsi.upc.edu
Email bix@sendu.me.uk

=cut

use strict;
use Bio::DB::Taxonomy;
use Bio::TreeIO;
use Bio::Tree::Compatible;
use Getopt::Long;

my @species;
my $index_dir = "./db/";
my $nodesfile = "nodes.dmp";
my $namesfile = "names.dmp";
my $use_entrez = 0;

# the input to the script is an array of species names
GetOptions('s|species=s' => \@species, 'd|dir:s' => \$index_dir, 'o|nodesfile:s' => \$nodesfile, 'a|namesfile:s' => \$namesfile, 'e|entrez' => \$use_entrez);


my $db = new Bio::DB::Taxonomy(-source => $use_entrez ? 'entrez' : 'flatfile',
                               -directory => $index_dir,
                               -nodesfile => $nodesfile,
                               -namesfile => $namesfile);

# the full lineages of the species are merged into a single tree
my $tree;
for my $name (@species) {
  my $ncbi_id = $db->get_taxonid($name);
  if ($ncbi_id) {
    my $node = $db->get_taxon(-taxonid => $ncbi_id);
    
    if ($tree) {
      $tree->merge_lineage($node);
    }
    else {
      $tree = new Bio::Tree::Tree(-node => $node);
    }
  }
  else {
    warn "no NCBI Taxonomy node for species ",$name,"\n";
  }
}

# simple paths are contracted by removing degree one nodes
$tree->contract_linear_paths;

# convert tree ids to their names for nice output with TreeIO
foreach my $node ($tree->get_nodes) {
  $node->id($node->node_name);
}

# the tree is output in Newick format
my $output = new Bio::TreeIO(-format => 'newick');
$output->write_tree($tree);

1;
