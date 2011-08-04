#!/usr/bin/perl
# This is a -*-Perl-* file (make my emacs happy)

=head1 NAME

bp_query_entrez_taxa - query Entrez taxonomy database and print out information 

=head1 USAGE

bp_query_entrez_taxa "Homo sapiens" "Saccharomyces cerevisiae" Rhizopus Metazoa
bp_query_entrez_taxa -gi 28800981 -gi 54301680 -db nucleotide
bp_query_entrez_taxa -gi 71836523 -db protein

 Provide the genus and species name in quotes, you can also query for
 a non-species node like Family or Order

Command-line options:
   -v or --verbose  : print verbose debugging info
   -gi              : one or many GI numbers to lookup taxon id for
   -db              : the sequence db (nucleotide or protein) the GI is for

   other arguments are assumed to be species names to lookup in taxonomy db


=head1 AUTHOR

Jason Stajich jason-at-bioperl-dot-org

=cut

use strict;
use warnings;
use Bio::DB::Taxonomy;
use Getopt::Long;

my $verbose = 0;
my (@gi, $dbname);
GetOptions('v|verbose' => \$verbose,
	   'gi:i'      => \@gi,
	   'db:s'      => \$dbname);

my $db = new Bio::DB::Taxonomy(-source => 'entrez', -verbose => $verbose);
if( @gi ) {
    my @nodes= $db->get_Taxonomy_Node(-gi => \@gi,
				      -db => $dbname);
    for my $node ( @nodes ) {
	my $gi = shift @gi;
	print " for gi $gi:\n";
	print " taxonid is ",$node->ncbi_taxid,"\n";    
	print " node is ", join(", ",$node->classification), "\n";
	print " species is ", $node->species,"\n";
	print " parent is ", $node->parent_id, "\n";
	print " rank is ", $node->rank, "\n";
	print " genetic_code  ", $node->genetic_code, "\n";
	print " mito_genetic_code  ", $node->mitochondrial_genetic_code, "\n";
	print " scientfic name is ", $node->binomial, "\n";
    }	
}

print "\n\n";
for my $name ( @ARGV ) {
    my $taxonid = $db->get_taxonid($name);
    my $node   = $db->get_Taxonomy_Node(-taxonid => $taxonid);
    print "taxonid is $taxonid\n";

    print " node is ", join(", ",$node->classification), "\n";
    print " species is ", $node->species,"\n";
    print " parent is ", $node->parent_id, "\n";
    print " rank is ", $node->rank, "\n";
    print " genetic_code  ", $node->genetic_code, "\n";
    print " mito_genetic_code  ", $node->mitochondrial_genetic_code, "\n";
    print " scientfic name is ", $node->binomial, "\n";
    print " common name is ", $node->common_name, "\n";
    print " create date is ", $node->create_date, "\n";
    print " update date is ", $node->update_date, "\n";
    print " pub date is ", ($node->pub_date || ''), "\n";
    print " variant is ", $node->variant, "\n";
    print " sub_species is ", $node->sub_species, "\n";
    print " organelle is ", $node->organelle, "\n";
    print " division is ", $node->division, "\n";
}
