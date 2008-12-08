# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 5);
    
	use_ok('Bio::Tree::RandomFactory');
    use_ok('Bio::TreeIO');
    use_ok('Bio::Tree::Statistics');
}

my $FILE1 = test_output_file();
 
my $ssize = 5;
my $factory = Bio::Tree::RandomFactory->new(-sample_size => $ssize);
my $stats = Bio::Tree::Statistics->new();

my $tree = $factory->next_tree;

is($tree->get_nodes, ($ssize * 2 - 1));

my $treeio = Bio::TreeIO->new(-format => 'newick', -file => ">$FILE1");

$treeio->write_tree($tree);
undef $treeio;

ok(-s $FILE1);
