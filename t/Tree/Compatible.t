# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
  use lib '.';
  use Bio::Root::Test;
  
  test_begin(-tests => 5,
			 -requires_module => 'Set::Scalar');
  
  use_ok('Bio::Tree::Compatible');
  use_ok('Bio::TreeIO');
}

# these tests are done with direct access to Bio::Tree::Compatible methods,
# instead of via creating a Bio::Tree::Compatible->new() object or similar...
# the docs seem to indicate that is normally possible? TODO?

my $in = Bio::TreeIO->new(-format => 'newick',
						  -fh     => \*DATA);

# the common labels of (((A,B)C,D),(E,F,G)); and ((A,B)H,E,(J,(K)G)I);
# are [A,B,E,G]

my $t1 = $in->next_tree;
my $t2 = $in->next_tree;
my $common = Bio::Tree::Compatible::common_labels($t1,$t2);
my $labels = Set::Scalar->new(qw(A B E G));
ok($common->is_equal($labels));

# the topological restrictions of (((A,B)C,D),(E,F,G)); and
# ((A,B)H,E,(J,(K)G)I); to their common labels, [A,B,E,G], are,
# respectively, ((A,B),(E,G)); and ((A,B),E,(G));

Bio::Tree::Compatible::topological_restriction($t1,$common);
Bio::Tree::Compatible::topological_restriction($t2,$common);
my $t3 = $in->next_tree;
my $t4 = $in->next_tree;
# ok($t1->is_equal($t3)); # is_equal method missing in Bio::Tree::Tree
# ok($t2->is_equal($t4)); # is_equal method missing in Bio::Tree::Tree

# the topological restrictions of (((A,B)C,D),(E,F,G)); and
# ((A,B)H,E,(J,(K)G)I); to their common labels, [A,B,E,G], are
# compatible

my ($incompat, $ilabels, $inodes) = Bio::Tree::Compatible::is_compatible($t3,$t4);
ok(!$incompat);

# (((B,A),C),D); and ((A,(D,B)),C); are incompatible

my $t5 = $in->next_tree;
my $t6 = $in->next_tree;
($incompat, $ilabels, $inodes) = Bio::Tree::Compatible::is_compatible($t5,$t6);
ok($incompat);

__DATA__
(((A,B)C,D),(E,F,G));
((A,B)H,E,(J,(K)G)I);
((A,B),(E,G));
((A,B),E,(G));
(((B,A),C),D);
((A,(D,B)),C);
