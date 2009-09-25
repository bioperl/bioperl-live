# -*-Perl-*- Test Harness script for Bioperl
# $Id: epost.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 65,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

# check -correspondence => 0 (default) 
my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'elink',
    -file       => test_input_file('eutils','elink_nhist.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Link');
is(join(',',$eutil->get_databases), 'pubmed');

# for elinks, IDs are globbed together when called from the parser unless a
# database is specified when cmd=neighbor_history is used, no IDs come back
# (they are stored on the server for further work)

is(join(',',$eutil->get_ids), '', 'get_ids');
my @ls = $eutil->get_LinkSets;
is(scalar(@ls), 2, 'uncorrelated LinkSets lump everything together');
is(join(',',$ls[1]->get_databases), 'pubmed');
isa_ok($ls[0], 'Bio::Tools::EUtilities::EUtilDataI');
isa_ok($ls[0], 'Bio::Tools::EUtilities::Link::LinkSet');
isa_ok($ls[0], 'Bio::Tools::EUtilities::HistoryI');

# check data in LinkSets
# Note that retrieved IDs and submitted IDs are lumped together (don't correspond)
is(join(',',$ls[0]->get_ids), '');
is(join(',',$ls[0]->get_databases), 'pubmed');
is(join(',',$ls[0]->get_submitted_ids), '730439,68536103,1621261,20807972');
is($ls[0]->get_dbfrom, 'protein');
is(join(',',$ls[0]->get_link_names), 'protein_pubmed');
is($ls[0]->has_scores, 0);
is($ls[0]->has_linkout, 0);
is($ls[0]->has_neighbor, 0);

# no LinkInfo
my @info = $ls[0]->get_LinkInfo;
is(scalar(@info), 0);

# no UrlLinks
my @urls = $ls[0]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[0]->get_webenv, '085LBC0s_G5ZenmRAnAm9dgF-TYrzyM9zVawz6_GfunjA5iasUqoGSfSzd@991070AE944054A1_0001SID');
is($ls[0]->get_query_key, 1);

# next 
is(join(',',$ls[1]->get_ids), '');
is(join(',',$ls[1]->get_databases), 'pubmed');
is(join(',',$ls[1]->get_submitted_ids), '730439,68536103,1621261,20807972');
is(join(',',$ls[1]->get_link_names), 'protein_pubmed_refseq');
is($ls[1]->get_dbfrom, 'protein');
is($ls[1]->has_scores, 0);
is($ls[1]->has_linkout, 0);
is($ls[1]->has_neighbor, 0);

# no LinkInfo
@info = $ls[1]->get_LinkInfo;
is(scalar(@info), 0);

# no UrlLinks
@urls = $ls[1]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[1]->get_webenv, '085LBC0s_G5ZenmRAnAm9dgF-TYrzyM9zVawz6_GfunjA5iasUqoGSfSzd@991070AE944054A1_0001SID');
is($ls[1]->get_query_key, 2);

# check -correspondence => 1
$eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'elink',
    -file       => test_input_file('eutils','elink_nhist_corr.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Link');
is(join(',',$eutil->get_databases), 'pubmed');

# for elinks, IDs are globbed together when called from the parser unless a database is specified
is(join(',',$eutil->get_ids), '', 'get_ids');
@ls = $eutil->get_LinkSets;
is(scalar(@ls), 6, 'correlated LinkSets separate ID data');
is(join(',',$ls[1]->get_databases), 'pubmed');
isa_ok($ls[0], 'Bio::Tools::EUtilities::EUtilDataI');
isa_ok($ls[0], 'Bio::Tools::EUtilities::Link::LinkSet');
isa_ok($ls[0], 'Bio::Tools::EUtilities::HistoryI');

# check data in LinkSets
# Note that you can get more that one returned ID, but only one submitted ID
is(join(',',$ls[0]->get_ids), ''); 
is(join(',',$ls[0]->get_submitted_ids), '1621261');
is(join(',',$ls[0]->get_link_names), 'protein_pubmed');
is($ls[0]->get_dbfrom, 'protein');
is($ls[0]->has_scores, 0);
is($ls[0]->has_linkout, 0);
is($ls[0]->has_neighbor, 0);

# no LinkInfo
@info = $ls[0]->get_LinkInfo;
is(scalar(@info), 0);

# no UrlLinks
@urls = $ls[0]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[0]->get_webenv, '0-g5Po62X-zBqwiLv9LDfH6dJvaMByxF-B7jUpwxS73UvKdcD2qdti4CNbY@03F16D1B94400731_0005SID');
is($ls[0]->get_query_key, 1);

is(join(',',$ls[1]->get_ids), '');
is(join(',',$ls[1]->get_databases), 'pubmed');
is(join(',',$ls[1]->get_submitted_ids), '68536103');
is($ls[1]->get_dbfrom, 'protein');
is(join(',',$ls[1]->get_link_names), 'protein_pubmed');
is($ls[1]->has_scores, 0);
is($ls[1]->has_linkout, 0);
is($ls[1]->has_neighbor, 0);

# no LinkInfo
@info = $ls[1]->get_LinkInfo;
is(scalar(@info), 0);

# no UrlLinks
@urls = $ls[1]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[1]->get_webenv, '0-g5Po62X-zBqwiLv9LDfH6dJvaMByxF-B7jUpwxS73UvKdcD2qdti4CNbY@03F16D1B94400731_0005SID');
is($ls[1]->get_query_key, 2);
