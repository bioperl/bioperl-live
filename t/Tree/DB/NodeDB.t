# -*-Perl-*- Test Harness script for Bioperl

use strict;

use Bio::Root::Test;
use Bio::DB::Tree::Node;

my $node = Bio::DB::Tree::Node->new(-node_id => 1);

is($node->node_id, 1);
done_testing();
exit;
