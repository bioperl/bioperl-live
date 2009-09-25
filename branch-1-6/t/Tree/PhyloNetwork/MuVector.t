# $Id$
use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;

  test_begin(-tests => 8);

  use_ok('Bio::PhyloNetwork::muVector');
}

my $vec1=Bio::PhyloNetwork::muVector->new(4);
my $vec2=Bio::PhyloNetwork::muVector->new([1,2,3,4]);
isa_ok($vec1,'Bio::PhyloNetwork::muVector');
isa_ok($vec1,'Bio::PhyloNetwork::muVector');

my $vec3=-1*$vec2;
my $vec4=$vec3+$vec2;

is($vec4 cmp $vec1,0,'arithmetic');
ok($vec2->display() eq "(1 2 3 4)",'display');
ok($vec2->is_positive(),'is_positive');

my $vec5=Bio::PhyloNetwork::muVector->new([2,3,5,0,77]);
my $vec6=Bio::PhyloNetwork::muVector->new([2,3,4,5,-7]);

ok($vec5->geq_poset($vec6) == 0,'geq_poset');

my $vec7=Bio::PhyloNetwork::muVector->new([2,3,5,0,77]);
my $vec8=Bio::PhyloNetwork::muVector->new([2,3,4,-1,-7]);

ok($vec7->geq_poset($vec8),'geq_poset');

