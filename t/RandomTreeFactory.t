# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

my $error;
use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    $error = 0; 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 3;
}

use Bio::Tree::RandomFactory;
use Bio::TreeIO;
ok(1);

#END { unlink("out.tre") }
 
my $ssize = 6;
my $factory = new Bio::Tree::RandomFactory(-sample_size => $ssize);

my $tree = $factory->next_tree;

ok($tree->get_nodes, ($ssize * 2 - 1));

my $treeio = new Bio::TreeIO(-format => 'newick', -file => ">out.tre");

$treeio->write_tree($tree);
ok(1);
