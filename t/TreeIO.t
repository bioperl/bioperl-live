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
    plan tests => 3; 

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

ok(1);

my $treeio = new Bio::TreeIO(-format => 'newick',
			     -file   => Bio::Root::IO->catfile('t','data', 
							       'LOAD_Ccd1.dnd'));
ok($treeio);
my $tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

my @nodes = $tree->get_nodes;

foreach my $node ( @nodes ) {
    if( $node->isa('Bio::Tree::PhyloNode') ) {
	print $node->id, " ", $node->bootstrap, "\n";
    } else { 

    }
}
# ok(@nodes, 5);
