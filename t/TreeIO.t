# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }

    use Test;
    plan tests => 4; 

#    eval { require XML::Parser::PerlSAX; };
#    if( $@ ) {
#	print STDERR "XML::Parser::PerlSAX not loaded. This means TreeIO::phyloxml test cannot be executed. Skipping\n";
#	foreach ( 1..43 ) {
#	    skip(1,1);
#	}
#       $error = 1;
#	
#    } 

}

if( $error == 1 ) {
    exit(0);
}

use Bio::TreeIO;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $treeio = new Bio::TreeIO(-verbose => $verbose,
			     -format => 'newick',
			     -file   => Bio::Root::IO->catfile('t','data', 
							       'cysprot1b.dnd'));
#							       'LOAD_Ccd1.dnd'));
ok($treeio);
my $tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

my @nodes = $tree->get_nodes;
ok(@nodes, 7);

foreach my $node ( @nodes ) {
    if( $node->isa('Bio::Tree::PhyloNode') ) {
#	if( $verbose ) { print "id=", $node->id, 
#			 "; branch_len=", $node->branch_length, "\n"; }
    } else { 
#	print "node was ", $node->to_string(), "\n" if( $verbose );
    }
}
$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format => 'newick', 
			  -file   => '>testnewick.phylip');
$treeio->write_tree($tree);

$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format => 'newick',
			  -file   => Bio::Root::IO->catfile('t','data', 
							    'LOAD_Ccd1.dnd'));
ok($treeio);
$tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

@nodes = $tree->get_nodes;
ok(@nodes, 53);

foreach my $node ( @nodes ) {
    if( $node->isa('Bio::Tree::PhyloNode') ) {
	print $node->id, ":", $node->branch_length(), "\n" if( $verbose );

#	if( $verbose ) { print "id=", $node->id, 
#			 "; branch_len=", $node->branch_length, "\n"; }
    } else { 
	print ":", $node->branch_length(), "\n" if( $verbose );
    }
}

