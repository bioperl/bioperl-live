# $Id$
# -*-Perl-*-

use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;

  test_begin(-tests => 70,
	     -requires_modules => [qw(Bio::PhyloNetwork
				      Bio::PhyloNetwork::TreeFactory)]);

  use_ok('Bio::PhyloNetwork::Factory');
}

my $factory = Bio::PhyloNetwork::Factory->new(-numleaves=>3);
isa_ok($factory,'Bio::PhyloNetwork::Factory');

my @nets;

while (my $net=$factory->next_network()) {
  push @nets,$net;
}

is(scalar @nets,66,'factory');

my @netsbk;

foreach my $enw (<DATA>) {
  my $net = Bio::PhyloNetwork->new(-eNewick=>$enw);
  push @netsbk,$net;
}
is(scalar @netsbk,66,'seen all the data lines');

foreach my $net (@nets) {
  my $found=0;
  foreach my $netbk (@netsbk) {
    if ($net->mu_distance($netbk)==0) {
      $found = 1;
      last;
    }
  }
  ok($found,'found');
}

__DATA__
((l1,l3),l2); 
(#H1,(#H1,(l1,l3))); (l2)#H1; 
(#H1,(#H1,((#H2,l1),#H2))); (l2)#H1; (l3)#H2; 
(#H1,(#H2,(#H1,(#H2,l1)))); (l2)#H1; (l3)#H2; 
(#H1,(#H1,(#H2,l1))); ((#H2,l2))#H1; (l3)#H2; 
(#H1,(#H1,((#H2,l3),#H2))); (l2)#H1; (l1)#H2; 
(#H1,(#H2,(#H1,(l3,#H2)))); (l2)#H1; (l1)#H2; 
(#H1,(#H1,(l3,#H2))); ((#H2,l2))#H1; (l1)#H2; 
(#H1,(l1,(l3,#H1))); (l2)#H1; 
((((l3,#H1),#H2),#H2),#H1); (l2)#H1; (l1)#H2; 
(#H1,((#H1,l3),#H2)); ((l2,#H2))#H1; (l1)#H2; 
(#H1,(#H2,(#H2,l1))); (l2)#H1; ((#H1,l3))#H2; 
(#H1,(l3,(l1,#H1))); (l2)#H1; 
((((l1,#H1),#H2),#H2),#H1); (l2)#H1; (l3)#H2; 
(#H1,((#H1,l1),#H2)); ((l2,#H2))#H1; (l3)#H2; 
(#H1,(#H2,(#H2,l3))); (l2)#H1; ((#H1,l1))#H2; 
(#H1,(l2,#H1)); ((l1,l3))#H1; 
((#H1,(l2,#H2)),#H1); ((#H2,l1))#H1; (l3)#H2; 
(#H1,(#H1,l2)); ((#H2,(#H2,l1)))#H1; (l3)#H2; 
(#H1,(#H2,(#H1,l2))); ((#H2,l1))#H1; (l3)#H2; 
((#H1,(l2,#H2)),#H1); ((#H2,l3))#H1; (l1)#H2; 
(#H1,(#H1,l2)); ((#H2,(#H2,l3)))#H1; (l1)#H2; 
(#H1,(#H2,(#H1,l2))); ((l3,#H2))#H1; (l1)#H2; 
((l2,#H1),(l1,#H1)); (l3)#H1; 
(l2,(#H1,(l1,#H1))); (l3)#H1; 
((l2,#H1),(l3,#H1)); (l1)#H1; 
(l2,(#H1,(l3,#H1))); (l1)#H1; 
(l1,(l3,l2)); 
(l1,(#H1,(#H1,l3))); (l2)#H1; 
((#H1,((#H1,l3),#H2)),#H2); (l2)#H1; (l1)#H2; 
((#H1,(#H1,(#H2,l3))),#H2); (l2)#H1; (l1)#H2; 
(#H2,(#H2,(#H1,(#H1,l3)))); (l2)#H1; (l1)#H2; 
((#H1,(l3,#H1)),#H2); ((#H2,l2))#H1; (l1)#H2; 
(#H2,(#H2,l1)); (l2)#H1; (((l3,#H1),#H1))#H2; 
((l1,#H1),(l3,#H1)); (l2)#H1; 
((#H1,(l1,#H2)),#H2); (l2)#H1; ((#H1,l3))#H2; 
(#H2,(#H2,(#H1,l1))); (l2)#H1; ((#H1,l3))#H2; 
(#H2,(#H2,(#H1,l3))); (l2)#H1; ((#H1,l1))#H2; 
(#H2,(#H1,(l3,#H2))); (l2)#H1; ((#H1,l1))#H2; 
(#H1,(l1,#H1)); ((l3,l2))#H1; 
(#H1,(l1,#H1)); ((#H2,(#H2,l2)))#H1; (l3)#H2; 
(l1,(#H1,(l2,#H1))); (l3)#H1; 
((#H1,(#H1,(#H2,l2))),#H2); (l3)#H1; (l1)#H2; 
(#H2,(((l2,#H1),#H1),#H2)); (l3)#H1; (l1)#H2; 
(#H2,(((l2,#H1),#H2),#H1)); (l3)#H1; (l1)#H2; 
((#H1,(#H1,l2)),#H2); ((l3,#H2))#H1; (l1)#H2; 
(#H1,(l3,(l2,#H1))); (l1)#H1; 
(#H1,(#H1,(l2,l3))); (l1)#H1; 
(#H1,((l3,#H1),l2)); (l1)#H1; 
((l1,l2),l3); 
(#H1,(#H1,(l1,l2))); (l3)#H1; 
(#H1,(#H1,((#H2,l1),#H2))); (l3)#H1; (l2)#H2; 
(#H1,(#H2,(#H1,(#H2,l1)))); (l3)#H1; (l2)#H2; 
(#H1,(#H1,((#H2,l2),#H2))); (l3)#H1; (l1)#H2; 
(#H1,(#H2,(#H1,(l2,#H2)))); (l3)#H1; (l1)#H2; 
(#H1,(l1,(l2,#H1))); (l3)#H1; 
((((l2,#H1),#H2),#H2),#H1); (l3)#H1; (l1)#H2; 
(#H1,(#H2,(#H2,l1))); (l3)#H1; ((#H1,l2))#H2; 
(#H1,(l2,(l1,#H1))); (l3)#H1; 
((((l1,#H1),#H2),#H2),#H1); (l3)#H1; (l2)#H2; 
(#H1,(#H2,(#H2,l2))); (l3)#H1; ((#H1,l1))#H2; 
(#H1,(l3,#H1)); ((l1,l2))#H1; 
(#H1,(#H1,l3)); ((#H2,(#H2,l1)))#H1; (l2)#H2; 
(#H1,(#H1,l3)); ((#H2,(#H2,l2)))#H1; (l1)#H2; 
(l3,(#H1,(l1,#H1))); (l2)#H1; 
(l3,(#H1,(l2,#H1))); (l1)#H1; 
