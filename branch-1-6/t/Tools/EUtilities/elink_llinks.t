# -*-Perl-*- Test Harness script for Bioperl
# $Id: elink_llinks.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 86,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

# check -correspondence => 0 (default) - this is set up to return the
# exact same thing as correspondece = 1, tested below)
my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'elink',
    -file       => test_input_file('eutils','elink_llinks.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Link');
is(join(',',$eutil->get_databases), '');

# for elinks, IDs are globbed together when called from the parser unless a database is specified
is(join(',',$eutil->get_ids), '730439,68536103,1621261,20807972', 'get_ids');
my @ls = $eutil->get_LinkSets;
is(scalar(@ls), 4, 'uncorrelated LinkSets lump everything together');
is(join(',',$ls[1]->get_databases), '');
isa_ok($ls[0], 'Bio::Tools::EUtilities::EUtilDataI');
isa_ok($ls[0], 'Bio::Tools::EUtilities::Link::LinkSet');

# check data in LinkSets
is(join(',',$ls[0]->get_ids), '730439');
is(join(',',$ls[0]->get_databases), '');
is(join(',',$ls[0]->get_submitted_ids), '730439');
is($ls[0]->get_dbfrom, 'protein');
is(join(',',$ls[0]->get_link_names), '');
is($ls[0]->has_scores, 0);
is($ls[0]->has_linkout, 0);
is($ls[0]->has_neighbor, 0);

# HistoryI
is($ls[0]->get_webenv, undef);
is($ls[0]->get_query_key, undef);

# no LinkInfo
my @info = $ls[0]->get_LinkInfo;
is(scalar(@info), 0);

# has UrlLink!
my @urls = $ls[0]->get_UrlLinks;
is(scalar(@urls), 1);

is($urls[0]->get_dbfrom, 'protein');
is($urls[0]->get_attribute, 'registration required');
is($urls[0]->get_icon_url, 'http://www.ncbi.nlm.nih.gov/entrez/query/egifs/http:--salilab.org-modbase-img-modbase.jpg');
is($urls[0]->get_subject_type, 'structure');
is($urls[0]->get_url, 'http://salilab.org/modbase/search?databaseID=730439');
is($urls[0]->get_link_name, undef);
is($urls[0]->get_provider_name, 'MODBASE, Database of Comparative Protein Structure Models (Sali Lab/UCSF)');
is($urls[0]->get_provider_abbr, 'modbase');
is($urls[0]->get_provider_id, 5447);
is($urls[0]->get_provider_url, 'http://salilab.org/modbase');
is($urls[0]->get_provider_icon_url, 'http://salilab.org/modbase/img/modbase.jpg');

# Another LinkSet...
is(join(',',$ls[1]->get_ids), '68536103');
is(join(',',$ls[1]->get_databases), '');
is(join(',',$ls[1]->get_submitted_ids), '68536103');
is(join(',',$ls[1]->get_link_names), '');
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
is($ls[1]->get_webenv, undef);
is($ls[1]->get_query_key, undef);

# check -correspondence => 1
$eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'elink',
    -file       => test_input_file('eutils','elink_llinks_corr.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Link');
is(join(',',$eutil->get_databases), '');

# for elinks, IDs are globbed together when called from the parser unless a database is specified
is(join(',',$eutil->get_ids), '1621261,68536103,20807972,730439', 'get_ids');
@ls = $eutil->get_LinkSets;
is(scalar(@ls), 4, 'correlated LinkSets separate ID data');
is(join(',',$ls[1]->get_databases), '');
isa_ok($ls[0], 'Bio::Tools::EUtilities::EUtilDataI');
isa_ok($ls[0], 'Bio::Tools::EUtilities::Link::LinkSet');

# check data in LinkSets
is(join(',',$ls[0]->get_ids), '1621261');
is(join(',',$ls[0]->get_databases), '');
is(join(',',$ls[0]->get_submitted_ids), '1621261');
is(join(',',$ls[0]->get_link_names), '');
is($ls[0]->get_dbfrom, 'protein');
is($ls[0]->has_scores, 0);
is($ls[0]->has_linkout, 0);
is($ls[0]->has_neighbor, 0);

# HistoryI
is($ls[0]->get_webenv, undef);
is($ls[0]->get_query_key, undef);

# no LinkInfo
@info = $ls[0]->get_LinkInfo;
is(scalar(@info), 0);

# has UrlLink!
@urls = $ls[0]->get_UrlLinks;
is(scalar(@urls), 1);

is($urls[0]->get_dbfrom, 'protein');
is($urls[0]->get_attribute, 'registration required');
is($urls[0]->get_icon_url, 'http://www.ncbi.nlm.nih.gov/entrez/query/egifs/http:--salilab.org-modbase-img-modbase.jpg');
is($urls[0]->get_subject_type, 'structure');
is($urls[0]->get_url, 'http://salilab.org/modbase/search?databaseID=1621261');
is($urls[0]->get_link_name, undef);
is($urls[0]->get_provider_name, 'MODBASE, Database of Comparative Protein Structure Models (Sali Lab/UCSF)');
is($urls[0]->get_provider_abbr, 'modbase');
is($urls[0]->get_provider_id, 5447);
is($urls[0]->get_provider_url, 'http://salilab.org/modbase');
is($urls[0]->get_provider_icon_url, 'http://salilab.org/modbase/img/modbase.jpg');

# another LinkSet
is(join(',',$ls[1]->get_ids), '68536103');
is(join(',',$ls[1]->get_databases), '');
is(join(',',$ls[1]->get_submitted_ids), '68536103');
is($ls[1]->get_dbfrom, 'protein');
is(join(',',$ls[1]->get_link_names), '');
is($ls[1]->has_scores, 0);
is($ls[1]->has_linkout, 0);
is($ls[1]->has_neighbor, 0);

# no LinkInfo
@info = $ls[1]->get_LinkInfo;
is(scalar(@info), 0);

# no UrlLink
@urls = $ls[1]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[1]->get_webenv, undef);
is($ls[1]->get_query_key, undef);