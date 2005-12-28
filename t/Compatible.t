# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

my $error;
use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

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
  use vars qw($TESTCOUNT);
  $TESTCOUNT = 5;
  plan tests => $TESTCOUNT;
}

use Bio::Tree::Compatible;
use Bio::TreeIO;
my $verbose = 0;

my $in = new Bio::TreeIO(-format => 'newick',
			 -fh     => \*DATA);

# the common labels of (((A,B)C,D),(E,F,G)); and ((A,B)H,E,(J,(K)G)I);
# are [A,B,E,G]

my $t1 = $in->next_tree;
my $t2 = $in->next_tree;
my $common = $t1->Bio::Tree::Compatible::common_labels($t2);
my $labels = Set::Scalar->new(qw(A B E G));
ok($common->is_equal($labels));

# the topological restrictions of (((A,B)C,D),(E,F,G)); and
# ((A,B)H,E,(J,(K)G)I); to their common labels, [A,B,E,G], are,
# respectively, ((A,B),(E,G)); and ((A,B),E,(G));

$t1->Bio::Tree::Compatible::topological_restriction($common);
$t2->Bio::Tree::Compatible::topological_restriction($common);
my $t3 = $in->next_tree;
my $t4 = $in->next_tree;
# ok($t1->is_equal($t3)); # is_equal method missing in Bio::Tree::Tree
# ok($t2->is_equal($t4)); # is_equal method missing in Bio::Tree::Tree

# the topological restrictions of (((A,B)C,D),(E,F,G)); and
# ((A,B)H,E,(J,(K)G)I); to their common labels, [A,B,E,G], are
# compatible

my ($incompat, $ilabels, $inodes) = $t3->Bio::Tree::Compatible::is_compatible($t4);
ok(!$incompat);

# (((B,A),C),D); and ((A,(D,B)),C); are incompatible

my $t5 = $in->next_tree;
my $t6 = $in->next_tree;
my ($incompat, $ilabels, $inodes) = $t5->Bio::Tree::Compatible::is_compatible($t6);
ok($incompat);

__DATA__
(((A,B)C,D),(E,F,G));
((A,B)H,E,(J,(K)G)I);
((A,B),(E,G));
((A,B),E,(G));
(((B,A),C),D);
((A,(D,B)),C);
