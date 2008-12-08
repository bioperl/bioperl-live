# $Id$
# -*-Perl-*-
use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;

  test_begin(-tests => 19,
	     -requires_modules => [qw(Bio::PhyloNetwork)]);

  use_ok('Bio::PhyloNetwork::TreeFactory');
}

my $factory=Bio::PhyloNetwork::TreeFactory->new(-numleaves=>4);
isa_ok($factory,'Bio::PhyloNetwork::TreeFactory');

my @nets;

while (my $net=$factory->next_network()) {
  push @nets,$net;
}

is(scalar @nets,15,'tree factory');

my @netsbk;

foreach my $enw (<DATA>) {
  my $net=Bio::PhyloNetwork->new(-eNewick=>$enw);
  push @netsbk,$net;
}
is(scalar @netsbk,15,'seen all the data lines');

foreach my $net (@nets) {
  my $found=0;
  foreach my $netbk (@netsbk) {
    if ($net->mu_distance($netbk)==0) {
      $found=1;
      last;
    }
  }
  ok($found,'found');
}

__DATA__
(l2,((l1,l4),l3)); 
((l1,l3),(l2,l4)); 
((l1,(l3,l4)),l2); 
(l4,(l2,(l1,l3))); 
((l4,(l3,l1)),l2); 
((l2,l3),(l1,l4)); 
(l1,((l2,l4),l3)); 
(l1,((l3,l4),l2)); 
(l4,(l1,(l2,l3))); 
(l1,((l3,l2),l4)); 
(l3,((l1,l4),l2)); 
((l1,(l2,l4)),l3); 
((l1,l2),(l3,l4)); 
((l4,(l2,l1)),l3); 
(l4,(l3,(l1,l2))); 
