# This is -*-Perl-*- code#
# Bioperl Test Harness Script for Modules#
# $Id: protgraph.t,v 1.1 2004/03/13 23:45:32 radams Exp
 use vars qw($NUMTESTS $DEBUG $ERROR);
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
}
exit 0 if $ERROR ==  1;
require Bio::Graph::ProteinGraph;
require Bio::Graph::IO;
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
ok  $node   = $g2->nodes_by_id('A64696');
ok $node->accession_number, 'A64696';
	

