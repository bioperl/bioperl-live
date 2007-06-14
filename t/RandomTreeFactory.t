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
use Bio::Tree::Statistics;
ok(1);

use vars qw($FILE1);
$FILE1 = 'out.tre';
END { unlink $FILE1; }
 
my $ssize = 5;
my $factory = Bio::Tree::RandomFactory->new(-sample_size => $ssize);
my $stats = Bio::Tree::Statistics->new();

my $tree = $factory->next_tree;

ok($tree->get_nodes, ($ssize * 2 - 1));

my $treeio = Bio::TreeIO->new(-format => 'newick', -file => ">$FILE1");

$treeio->write_tree($tree);
undef $treeio;

ok(-s $FILE1);
