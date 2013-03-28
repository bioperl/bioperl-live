# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 2);

    use_ok('Bio::Cluster::UniGene');
}

my $clu = Bio::Cluster::UniGene->new();
isa_ok($clu, "Bio::AnnotatableI");
