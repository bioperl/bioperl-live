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
    plan tests => 4;
}

use Bio::Tree::RandomFactory;
use Bio::TreeIO;
use Bio::Tree::Statistics;
ok(1);

use vars qw($FILE1);
$FILE1 = 'out.tre';
END { unlink $FILE1; }
 
my $ssize = 10;
my $factory = new Bio::Tree::RandomFactory(-sample_size => $ssize);

my $tree = $factory->next_tree;

ok($tree->get_nodes, ($ssize * 2 - 1));

my $treeio = new Bio::TreeIO(-format => 'newick', -file => ">$FILE1");

$treeio->write_tree($tree);
undef $treeio;

ok(-s $FILE1);
my $mutcount = 100;
$factory->add_Mutations($tree, $mutcount);

my $stats = new Bio::Tree::Statistics();
my $D = $stats->fu_and_li_D($tree);
ok($D);

