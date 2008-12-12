# -*-Perl-*- Test Harness script for Bioperl
# $Id: egquery.t 15112 2008-12-08 18:12:38Z sendu $
#

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 19);
	
    use_ok('Bio::Tools::EUtilities');
}

# Normal esearch
my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'egquery',
    -file       => test_input_file('eutils','egquery.xml'));

is($eutil->get_db, 'pubmed', 'get_db');
is($eutil->get_database, 'pubmed', 'get_database');
is(scalar($eutil->get_databases), 35, 'get_databases');
is($eutil->get_term, 'Notch AND Mus musculus','get_term');

## eveything else undef or 0
is ($eutil->get_count('pubmed'), 1803, 'get_count');
is ($eutil->get_count('protein'), 534, 'get_count');
is ($eutil->get_count('cdd'), 0, 'get_count');

my @qs = $eutil->get_GlobalQueries;
is(scalar(@qs), 35, 'get_GlobalQueries');
is($qs[2]->get_term, 'Notch AND Mus musculus', 'get_term');
is($qs[2]->get_database, 'journals', 'get_term');
is($qs[2]->get_count, 0, 'get_term');
is($qs[2]->get_status, 'Term or Database is not found', 'get_term');
is($qs[2]->get_menu_name, 'Journals', 'get_term');

is($qs[20]->get_term, 'Notch AND Mus musculus', 'get_term');
is($qs[20]->get_database, 'unists', 'get_term');
is($qs[20]->get_count, 61, 'get_term');
is($qs[20]->get_status, 'Ok', 'get_term');
is($qs[20]->get_menu_name, 'UniSTS', 'get_term');
