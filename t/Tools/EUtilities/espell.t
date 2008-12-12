# -*-Perl-*- Test Harness script for Bioperl
# $Id: espell.t 15112 2008-12-08 18:12:38Z sendu $
#

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 22,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

# Normal esearch
my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'espell',
    -file       => test_input_file('eutils','espell.xml'));

is($eutil->get_db, 'pubmed', 'get_db');
is(($eutil->get_dbs)[0], 'pubmed', 'get_dbs');
is($eutil->get_database, 'pubmed', 'get_database');
is(($eutil->get_databases)[0], 'pubmed', 'get_databases');
is($eutil->get_term, 'Netch AND Mus musclus','get_term');
is($eutil->get_corrected_query, 'notch AND mus musculus' ,'get_corrected_query');
is(scalar($eutil->get_replaced_terms), 2,'get_replaced_terms');
is(join(',',$eutil->get_replaced_terms), 'notch,musculus','get_replaced_terms');

# eveything else undef or 0
is ($eutil->get_count, undef, 'get_count');
my $history = $eutil->next_History;
is($history, undef);
my @ids2 = $eutil->get_ids;
is(scalar(@ids2), 0, 'get_ids');
is($eutil->get_retstart, undef,'get_retstart');
is($eutil->get_retmax, undef,'get_retmax');
is($eutil->get_translation_from, undef,'get_translation_from');
is($eutil->get_translation_to, undef,'get_translation_to');

# add Parameters
my $pb = Bio::Tools::EUtilities::EUtilParameters->new(-eutil => 'espell',
									   -db => 'protein',
									   -term => 'Notch AND Mus musculus');

is($eutil->get_db, 'pubmed', 'get_db');
is(($eutil->get_dbs)[0], 'pubmed', 'get_dbs');
is($eutil->get_database, 'pubmed', 'get_database');
is(($eutil->get_databases)[0], 'pubmed', 'get_databases');
is($eutil->get_term, 'Netch AND Mus musclus','get_term');