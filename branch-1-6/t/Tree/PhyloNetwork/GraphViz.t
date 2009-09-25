# $Id$
use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;

  test_begin(-tests => 3,
	     -requires_modules => [qw(Bio::PhyloNetwork
				      GraphViz)]);

  use_ok('Bio::PhyloNetwork::GraphViz');
}

my $net=Bio::PhyloNetwork->new(-eNewick=>'((H1,(H1,(H2,l))),H2)t0; (some long label)H1; ("quoted label")H2;');

my $gv=Bio::PhyloNetwork::GraphViz->new(-net=>$net,-short_labels=>1);
isa_ok($gv,'Bio::PhyloNetwork::GraphViz');
isa_ok($gv,'GraphViz');
