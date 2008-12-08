# -*-Perl-*- Test Harness script for Bioperl
# $Id: RandomTreeFactory.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 13);

    use_ok('Bio::TreeIO');
    use_ok('Bio::Tree::Statistics');
}


use Data::Dumper;
my $in = Bio::TreeIO->new(-format => 'nexus',
			  -file   => test_input_file('traittree.nexus'));
my $tree = $in->next_tree;
my $node = $tree->find_node(-id => 'N14');


my $stats = Bio::Tree::Statistics->new();
is $stats->cherries($tree), 8, 'cherries';
is $stats->cherries($tree, $node), 4, 'cherries';

# traits
my $key = $tree->add_trait(test_input_file('traits.tab'), 3);
is ($key, 'intermediate', 'read traits');

is $stats->ps($tree, $key), 4, 'parsimony score';
is $stats->ps($tree, $key, $node), 1, 'subtree parsimony score';

is $stats->ai($tree, $key), 0.628906, 'association index';
is $stats->ai($tree, $key, $node), 0.062500, 'subtree association index';

my $mc = $stats->mc($tree, $key);
is ($mc->{blue}, 2, 'monophyletic clade size');
is ($mc->{red}, 4, 'monophyletic clade size');
$node = $tree->find_node(-id => 'N10');
$mc = $stats->mc($tree, $key, $node);
is ($mc->{blue}, 2, 'monophyletic clade size');
is ($mc->{red}, 2, 'monophyletic clade size');
