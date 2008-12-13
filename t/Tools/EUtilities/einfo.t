# -*-Perl-*- Test Harness script for Bioperl
# $Id: einfo.t 15112 2008-12-08 18:12:38Z sendu $
#

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 51,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

## einfo (no dbs)
my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'einfo',
    -file       => test_input_file('eutils','einfo_dbs.xml'));

is(scalar($eutil->get_available_databases), 37, 'get_available_databases');
is(scalar($eutil->get_databases), 37, 'get_databases');
is($eutil->get_db, 'pubmed', 'get_db');

# no data present for these
is($eutil->get_record_count, undef, 'get_record_count');
is($eutil->get_menu_name, undef, 'get_menu_name');
is($eutil->get_last_update, undef, 'get_last_update');
is($eutil->get_description, undef, 'get_description');

my @fields = $eutil->get_FieldInfo;
is(scalar(@fields), 0, 'FieldInfo');
my @linkinfo = $eutil->get_LinkInfo;
is(scalar(@linkinfo), 0, 'LinkInfo');

# einfo (db-specific)
$eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'einfo',
    -file       => test_input_file('eutils','einfo.xml'));

is($eutil->get_db, 'pubmed', 'get_db');
is(scalar($eutil->get_dbs), 1, 'get_dbs');
is($eutil->get_record_count, 18525976, 'get_record_count');
is($eutil->get_menu_name, 'PubMed', 'get_menu_name');
is($eutil->get_last_update, '2008/12/11 02:33', 'get_last_update');
is($eutil->get_description, 'PubMed bibliographic record', 'get_description');

@fields = $eutil->get_FieldInfo;
is(scalar(@fields), 41, 'FieldInfo');
# test two
is($fields[1]->get_term_count, 0, 'get_term_count');
is($fields[1]->get_field_name, 'UID', 'get_field_name');
is($fields[1]->get_field_code, 'UID', 'get_field_code');
is($fields[1]->get_field_description, 'Unique number assigned to publication', 'get_field_description');
is($fields[1]->is_date, 0, 'is_date');
is($fields[1]->is_singletoken, 1, 'is_singletoken');
is($fields[1]->is_hierarchy, 0, 'is_hierarchy');
is($fields[1]->is_hidden, 1, 'is_hidden');
is($fields[1]->is_numerical, 1, 'is_numerical');

is($fields[19]->get_term_count, 83, 'get_term_count');
is($fields[19]->get_field_name, 'MeSH Subheading', 'get_field_name');
is($fields[19]->get_field_code, 'SUBH', 'get_field_code');
is($fields[19]->get_field_description, 'Additional specificity for MeSH term', 'get_field_description');
is($fields[19]->is_date, 0, 'is_date');
is($fields[19]->is_singletoken, 1, 'is_singletoken');
is($fields[19]->is_hierarchy, 0, 'is_hierarchy');
is($fields[19]->is_hidden, 0, 'is_hidden');
is($fields[19]->is_numerical, 0, 'is_numerical');

@linkinfo = $eutil->get_LinkInfo;
is(scalar(@linkinfo), 46, 'LinkInfo');
# test two
is($linkinfo[1]->get_dbto, 'cancerchromosomes', 'get_dbto');
is($linkinfo[1]->get_dbfrom, 'pubmed', 'get_dbfrom');
is($linkinfo[1]->get_link_name, 'pubmed_cancerchromosomes', 'get_link_name');
is($linkinfo[1]->get_link_description, 'Related Cancer Chromosomes', 'get_link_description');
is($linkinfo[1]->get_priority, undef, 'get_priority');
is($linkinfo[1]->get_html_tag, undef, 'get_html_tag');
is($linkinfo[1]->get_url, undef, 'get_url');

is($linkinfo[12]->get_dbto, 'geo', 'get_dbto');
is($linkinfo[12]->get_dbfrom, 'pubmed', 'get_dbfrom');
is($linkinfo[12]->get_link_name, 'pubmed_geo', 'get_link_name');
is($linkinfo[12]->get_link_description, 'GEO records associated with pubmed record', 'get_link_description');
is($linkinfo[12]->get_priority, undef, 'get_priority');
is($linkinfo[12]->get_html_tag, undef, 'get_html_tag');
is($linkinfo[12]->get_url, undef, 'get_url');

