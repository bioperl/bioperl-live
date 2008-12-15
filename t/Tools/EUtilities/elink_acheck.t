# -*-Perl-*- Test Harness script for Bioperl
# $Id: epost.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 130,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

# check -correspondence => 0 (default) - this is set up to return the
# exact same thing as correspondece = 1, tested below)
my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'elink',
    -file       => test_input_file('eutils','elink_acheck.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Link');
is(join(',',$eutil->get_databases), 'LinkOut,cdd,gene,genome,genomeprj,'.
   'nuccore,pmc,protein,proteinclusters,pubmed,structure,taxonomy');

# for elinks, IDs are globbed together when called from the parser unless a database is specified
is(join(',',$eutil->get_ids('cdd')), '730439,68536103,1621261,20807972', 'get_ids');
is(join(',',$eutil->get_ids('LinkOut')), '730439,1621261,20807972', 'get_ids');
my @ls = $eutil->get_LinkSets;
is(scalar(@ls), 4, 'uncorrelated LinkSets lump everything together');
is(join(',',$ls[1]->get_databases), 'cdd,gene,genome,genomeprj,nuccore,pmc,'.
   'protein,proteinclusters,pubmed,structure,taxonomy');
isa_ok($ls[0], 'Bio::Tools::EUtilities::EUtilDataI');
isa_ok($ls[0], 'Bio::Tools::EUtilities::Link::LinkSet');

# check data in LinkSets
is(join(',',$ls[0]->get_ids), '730439');
is(join(',',$ls[0]->get_databases), 'LinkOut,cdd,pmc,protein,pubmed,structure,'.
   'taxonomy');
is(join(',',$ls[0]->get_submitted_ids), '730439');
is($ls[0]->get_dbfrom, 'protein');
is(join(',',$ls[0]->get_link_names), 'protein_cdd,protein_cdd_concise_2,'.
   'protein_cdd_summary,protein_pmc,protein_protein,'.
   'protein_protein_cdart_summary,protein_protein_identical,protein_pubmed,'.
   'protein_structure_related,protein_taxonomy,ExternalLink');
is($ls[0]->has_scores, 0);
is($ls[0]->has_linkout, 1);
is($ls[0]->has_neighbor, 0);
my @info = $ls[0]->get_LinkInfo;
is(scalar(@info), 11);
is($info[1]->get_database, 'cdd');
is($info[1]->get_dbfrom, 'protein');
is($info[1]->get_link_name, 'protein_cdd_concise_2');
is($info[1]->get_link_description, undef);
is($info[1]->get_link_menu_name, 'Concise Conserved Domain Links');
is($info[1]->get_priority, 128);
is($info[1]->get_html_tag, undef);
is($info[1]->get_url, undef);

is($info[10]->get_database, 'LinkOut');
is($info[10]->get_dbfrom, 'protein');
is($info[10]->get_link_name, 'ExternalLink');
is($info[10]->get_link_description, undef);
is($info[10]->get_link_menu_name, 'LinkOut');
is($info[10]->get_priority, 255);
is($info[10]->get_html_tag, 'LinkOut');
is($info[10]->get_url, undef);

# no UrlLinks
my @urls = $ls[0]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[0]->get_webenv, undef);
is($ls[0]->get_query_key, undef);

is(join(',',$ls[1]->get_ids), '68536103');
is(join(',',$ls[1]->get_databases), 'cdd,gene,genome,genomeprj,nuccore,pmc,'.
   'protein,proteinclusters,pubmed,structure,taxonomy');
is(join(',',$ls[1]->get_submitted_ids), '68536103');
is(join(',',$ls[1]->get_link_names), 'protein_cdd,protein_cdd_concise_2,'.
   'protein_cdd_summary,protein_gene,protein_genome,protein_genomeprj,'.
   'protein_nuccore,protein_pmc,protein_protein,protein_protein_cdart_summary,'.
   'protein_protein_identical,protein_proteinclusters,protein_pubmed,'.
   'protein_pubmed_refseq,protein_structure_related,protein_taxonomy');
is($ls[1]->get_dbfrom, 'protein');
is($ls[1]->has_scores, 0);
is($ls[1]->has_linkout, 0);
is($ls[1]->has_neighbor, 0);
@info = $ls[1]->get_LinkInfo;
is(scalar(@info), 16);
is($info[1]->get_database, 'cdd');
is($info[1]->get_dbfrom, 'protein');
is($info[1]->get_link_name, 'protein_cdd_concise_2');
is($info[1]->get_link_description, undef);
is($info[1]->get_link_menu_name, 'Concise Conserved Domain Links');
is($info[1]->get_priority, 128);
is($info[1]->get_html_tag, undef);
is($info[1]->get_url, undef);
is($info[14]->get_database, 'structure');
is($info[14]->get_dbfrom, 'protein');
is($info[14]->get_link_name, 'protein_structure_related');
is($info[14]->get_link_description, undef);
is($info[14]->get_link_menu_name, undef);
is($info[14]->get_priority, 128);
is($info[14]->get_html_tag, 'Related Structure');
# Note the UID tag at end
is($info[14]->get_url, 'http://structure.ncbi.nlm.nih.gov/Structure/cblast/'.
   'cblast.cgi?client=entrez&query_gi=<@UID@>');

# no UrlLinks
@urls = $ls[1]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[1]->get_webenv, undef);
is($ls[1]->get_query_key, undef);

# check -correspondence => 1
$eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'elink',
    -file       => test_input_file('eutils','elink_acheck.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Link');
is(join(',',$eutil->get_databases), 'LinkOut,cdd,gene,genome,genomeprj,'.
   'nuccore,pmc,protein,proteinclusters,pubmed,structure,taxonomy');

# for elinks, IDs are globbed together when called from the parser unless a database is specified
is(join(',',$eutil->get_ids('cdd')), '730439,68536103,1621261,20807972', 'get_ids');
is(join(',',$eutil->get_ids('LinkOut')), '730439,1621261,20807972', 'get_ids');
@ls = $eutil->get_LinkSets;
is(scalar(@ls), 4, 'correlated LinkSets separate ID data');
is(join(',',$ls[1]->get_databases), 'cdd,gene,genome,genomeprj,nuccore,pmc,'.
   'protein,proteinclusters,pubmed,structure,taxonomy');
isa_ok($ls[0], 'Bio::Tools::EUtilities::EUtilDataI');
isa_ok($ls[0], 'Bio::Tools::EUtilities::Link::LinkSet');

# check data in LinkSets
is(join(',',$ls[0]->get_ids), '730439');
is(join(',',$ls[0]->get_databases), 'LinkOut,cdd,pmc,protein,pubmed,structure,'.
   'taxonomy');
is(join(',',$ls[0]->get_submitted_ids), '730439');
is(join(',',$ls[0]->get_link_names), 'protein_cdd,protein_cdd_concise_2,'.
   'protein_cdd_summary,protein_pmc,protein_protein,'.
   'protein_protein_cdart_summary,protein_protein_identical,'.
   'protein_pubmed,protein_structure_related,protein_taxonomy,ExternalLink');
is($ls[0]->get_dbfrom, 'protein');
is($ls[0]->has_scores, 0);
is($ls[0]->has_linkout, 1);
is($ls[0]->has_neighbor, 0);
@info = $ls[0]->get_LinkInfo;
is(scalar(@info), 11);
is($info[1]->get_database, 'cdd');
is($info[1]->get_dbfrom, 'protein');
is($info[1]->get_link_name, 'protein_cdd_concise_2');
is($info[1]->get_link_description, undef);
is($info[1]->get_link_menu_name, 'Concise Conserved Domain Links');
is($info[1]->get_priority, 128);
is($info[1]->get_html_tag, undef);
is($info[1]->get_url, undef);

is($info[10]->get_database, 'LinkOut');
is($info[10]->get_dbfrom, 'protein');
is($info[10]->get_link_name, 'ExternalLink');
is($info[10]->get_link_description, undef);
is($info[10]->get_link_menu_name, 'LinkOut');
is($info[10]->get_priority, 255);
is($info[10]->get_html_tag, 'LinkOut');
is($info[10]->get_url, undef);

# no UrlLinks
@urls = $ls[0]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[0]->get_webenv, undef);
is($ls[0]->get_query_key, undef);

is(join(',',$ls[1]->get_ids), '68536103');
is(join(',',$ls[1]->get_databases), 'cdd,gene,genome,genomeprj,nuccore,pmc,'.
   'protein,proteinclusters,pubmed,structure,taxonomy');
is(join(',',$ls[1]->get_submitted_ids), '68536103');
is($ls[1]->get_dbfrom, 'protein');
is(join(',',$ls[1]->get_link_names), 'protein_cdd,protein_cdd_concise_2,'.
   'protein_cdd_summary,protein_gene,protein_genome,protein_genomeprj,'.
   'protein_nuccore,protein_pmc,protein_protein,protein_protein_cdart_summary,'.
   'protein_protein_identical,protein_proteinclusters,protein_pubmed,'.
   'protein_pubmed_refseq,protein_structure_related,protein_taxonomy');
is($ls[1]->has_scores, 0);
is($ls[1]->has_linkout, 0);
is($ls[1]->has_neighbor, 0);
@info = $ls[1]->get_LinkInfo;
is(scalar(@info), 16);
is($info[1]->get_database, 'cdd');
is($info[1]->get_dbfrom, 'protein');
is($info[1]->get_link_name, 'protein_cdd_concise_2');
is($info[1]->get_link_description, undef);
is($info[1]->get_link_menu_name, 'Concise Conserved Domain Links');
is($info[1]->get_priority, 128);
is($info[1]->get_html_tag, undef);
is($info[1]->get_url, undef);
is($info[14]->get_database, 'structure');
is($info[14]->get_dbfrom, 'protein');
is($info[14]->get_link_name, 'protein_structure_related');
is($info[14]->get_link_description, undef);
is($info[14]->get_link_menu_name, undef);
is($info[14]->get_priority, 128);
is($info[14]->get_html_tag, 'Related Structure');
# Note the UID tag at end
is($info[14]->get_url, 'http://structure.ncbi.nlm.nih.gov/Structure/cblast/'.
   'cblast.cgi?client=entrez&query_gi=<@UID@>');

# no UrlLinks
@urls = $ls[1]->get_UrlLinks;
is(scalar(@urls), 0);

# HistoryI
is($ls[1]->get_webenv, undef);
is($ls[1]->get_query_key, undef);

