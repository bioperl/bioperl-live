# This is -*-Perl-*- code
# $Id$
use strict;

BEGIN { 
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    plan tests => 8;
}

use Bio::DB::Taxonomy;

my $db = new Bio::DB::Taxonomy(-source => 'entrez');

ok($db);


my ($taxaid) = $db->get_taxaid('homo sapiens');
ok($taxaid, 9606);

my $n = $db->get_Taxonomy_Node($taxaid);
ok($n);
unless( $n ) {
    for ( 1..4 ) { skip(1,'no species object could be created') }
} else {
    ok($n->species, 'sapiens');
    ok($n->common_name, 'human');
    ok($n->ncbi_taxid, 9606);
    ok($n->division, 'mammals');
}
my $yeastid = $db->get_taxaid('Saccharomyces cerevisiae');
ok($yeastid, 4932);

