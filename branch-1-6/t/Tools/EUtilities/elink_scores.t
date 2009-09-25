# -*-Perl-*- Test Harness script for Bioperl
# $Id: epost.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 58,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

# check -correspondence => 0 (default) - this is set up to return the
# exact same thing as correspondece = 1, tested below)
my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'elink',
    -file       => test_input_file('eutils','elink_scores.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Link');
is(join(',',$eutil->get_databases), 'protein');

# for elinks, IDs are globbed together when called from the parser unless a database is specified
is(join(',',$eutil->get_ids), '15622530,15921743,70607303,68567951,145702933,'.
   '146304683,6015889,13813749,15897502,15622530,74573864,15921743', 'get_ids');
my @ls = $eutil->get_LinkSets;
is(scalar(@ls), 2, 'uncorrelated LinkSets lump everything together');
is(join(',',$ls[1]->get_databases), 'protein');
isa_ok($ls[0], 'Bio::Tools::EUtilities::EUtilDataI');
isa_ok($ls[0], 'Bio::Tools::EUtilities::Link::LinkSet');

# check data in LinkSets
is(join(',',$ls[0]->get_ids), '15622530,15921743,70607303,68567951,145702933,'.
   '146304683,6015889,13813749,15897502');
is(join(',',$ls[0]->get_databases), 'protein');
is(join(',',$ls[0]->get_submitted_ids), '15622530');
is($ls[0]->get_dbfrom, 'protein');
is(join(',',$ls[0]->get_link_names), 'protein_protein');
is($ls[0]->has_linkout, 0);
is($ls[0]->has_neighbor, 0);

# has relatedness scores!
is($ls[0]->has_scores, 1);

my %sd = (
    15622530 => 2147483647,
    15921743 => 381,
    70607303 => 178,
    68567951 => 178,
    145702933 => 161,
    146304683 => 161,
    6015889 => 142,
    13813749 => 142,
    15897502 => 142);

my %sc = $ls[0]->get_scores;
for my $id ($ls[0]->get_ids) {
    ok(exists($sc{$id}));
    is($sc{$id}, $sd{$id});
    delete $sd{$id};
}
is(keys %sd, 0);

# no LinkInfo
my @info = $ls[0]->get_LinkInfo;
is(scalar(@info), 0);

# no UrlLinks
my @urls = $ls[0]->get_UrlLinks;
is(scalar(@urls), 0);

is(join(',',$ls[1]->get_ids), '15622530,74573864,15921743');
is(join(',',$ls[1]->get_databases), 'protein');
is(join(',',$ls[1]->get_submitted_ids), '15622530');
is(join(',',$ls[1]->get_link_names), 'protein_protein_identical');
is($ls[1]->get_webenv, undef);
is($ls[1]->get_dbfrom, 'protein');
is($ls[1]->has_linkout, 0);
is($ls[1]->has_neighbor, 0);

# has relatedness scores!
is($ls[1]->has_scores, 1);

%sd = (
    15622530 => 2147483647,
    74573864 => 0,
    15921743 => 0,
);

%sc = $ls[1]->get_scores;
for my $id ($ls[1]->get_ids) {
    ok(exists($sc{$id}));
    is($sc{$id}, $sd{$id});
    delete $sd{$id};
}

is(keys %sd, 0);

# HistoryI
is($ls[1]->get_webenv, undef);
is($ls[1]->get_query_key, undef);

# no LinkInfo
@info = $ls[1]->get_LinkInfo;
is(scalar(@info), 0);

# no UrlLinks
@urls = $ls[1]->get_UrlLinks;
is(scalar(@urls), 0);

