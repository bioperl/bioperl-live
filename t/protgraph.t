# This is -*-Perl-*- code#
# Bioperl Test Harness Script for Modules#
# $Id: protgraph.t,v 1.1 2004/03/13 23:45:32 radams Exp
 use vars qw($NUMTESTS $DEBUG $ERROR $XML_ERROR);
 use strict;
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {  
  # to handle systems with no installed Test module 
   # we include the t dir (where a copy of Test.pm is located)   
 # as a fallback  
  eval { require Test;}; 
   $ERROR = 0;   
    if( $@ ) {
	 use lib 't'; 
    }   
   use Test;  
  $NUMTESTS = 14;   
   plan tests => $NUMTESTS;
    eval {	require Class::AutoClass;	
         	require Clone;    
      }; 
   if( $@ ) {     
	   warn("Class::AutoClass or Clone not installed. " .
            " This means that the module is not usable. Skipping tests");
		$ERROR = 1;  
       }
	
    eval {
		require XML::Twig;
		};
	if ($@) {
		warn "XML::Twig needed for XML format parsing, skipping these tests";
		$XML_ERROR = 1;
		}
}
exit 0 if $ERROR ==  1;
require Bio::Graph::ProteinGraph;
require Bio::Graph::IO;
require Bio::Graph::Edge;
my $verbose = 0;
$verbose    = 1 if $DEBUG;
ok 1;
################1st of all let's test the io.....

###############  test dip tab delimited  format  ###################
## test read...
ok my $io = Bio::Graph::IO->new(-format => 'dip',
                                -file   => Bio::Root::IO->catfile("t","data","tab1part.mif"));
ok my $gr = $io->next_network();
ok my $node   = $gr->nodes_by_id('A64696');
ok $node->accession_number, 'A64696';

##test write. to filehandle...

ok my $out =  Bio::Graph::IO->new(-format => 'dip',
                                  -file   =>">". Bio::Root::IO->catfile("t","data","out.mif"));
ok $out->write_network($gr);  

## can we round trip, is out format same as original format?
ok my $io2 = Bio::Graph::IO->new(-format    => 'dip',
                                 -file     => Bio::Root::IO->catfile("t","data","out.mif"));
ok	my $g2     = $io2->next_network(); 
ok  $node      = $g2->nodes_by_id('A64696');
ok $node->accession_number, 'A64696';

############## now lets test the XML format.....
if (!$XML_ERROR){


}


##### now lets test some graph properties.....##
##before deleting there are 2
my @n1 = $gr->neighbors($gr->nodes_by_id('3047N'));
ok scalar @n1,2; 

ok $gr->remove_nodes($gr->nodes_by_id('3048N'));

## after deleting there is only 1 interactor
 @n1 = $gr->neighbors($gr->nodes_by_id('3047N'));
ok scalar @n1,1; 

##check no undefs left after node removal ##

ok map {$_->object_id}$gr->edges;
ok map {$_->object_id}$gr->nodes;





## count neighbors
my @ns = $gr->neighbors($gr->nodes_by_id('H64521'));
ok scalar @ns, 13;

##get CC using protein object..
ok  sprintf("%.3f", $gr->clustering_coefficient($gr->nodes_by_id('B64525'))), 0.067;

#.. andusing id string (same as previous, for convenience	)
ok  sprintf("%.3f", $gr->clustering_coefficient('B64525')), 0.067;	

## test has_node() method
ok $gr->has_node('B64525'), 1;
ok $gr->has_node('B64'), 0;

ok $gr->remove_dup_edges($gr->nodes_by_id('3103N'));
						


## remove dups
ok $gr->remove_dup_edges();
	
## should now be no duplicates
my @dupids = map{$_->object_id()} $gr->dup_edges();
ok $dupids[0], undef;

