# -*-Perl-*- Test Harness script for Bioperl
# $Id: esearch.t 15112 2008-12-08 18:12:38Z sendu $
#

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 33,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

my @ids = qw(6679096 31543332 134288853 483581 20805941 187951953 169158074
123228044 148676374 114326469 148707003 187952787 123233807 148694865 148694864
148694863 148694861 148694862 8705244 8568086);

# test any Query-related methods (term related)

# Normal esearch
my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'esearch',
    -file       => test_input_file('eutils','esearch1.xml'));

# w/o a ParameterBase, only IDs, count, retstart/retmax, optionally History
is ($eutil->get_count, 534, 'get_count');
my $history = $eutil->next_History;
is($history, undef);
my @ids2 = $eutil->get_ids;
is_deeply(\@ids2, \@ids, 'get_ids');
is($eutil->get_retstart, 0,'get_retstart');
is($eutil->get_retmax, 20,'get_retmax');
is($eutil->get_translation_from, 'Mus musculus','get_translation_from');
is($eutil->get_translation_to, '("Mus musculus"[Organism:__txid10090] OR Mus musculus)','get_translation_to');

# the database isn't carried into the parsers unless a EUtilParameters is present
is($eutil->get_db, undef, 'get_db');
is($eutil->get_database, undef, 'get_database');
is($eutil->get_term, undef,'get_term');

# add Parameters
my $pb = Bio::Tools::EUtilities::EUtilParameters->new(-eutil => 'esearch',
									   -db => 'protein',
									   -term => 'Notch AND Mus musculus');

$eutil->parameter_base($pb);

# now will work...
is($eutil->get_db, 'protein', 'get_db');
is($eutil->get_database, 'protein', 'get_database');
is($eutil->get_term, 'Notch AND Mus musculus','get_term');

# espell only (should be undef)
is($eutil->get_corrected_query, undef ,'get_corrected_query');
is($eutil->get_replaced_terms, undef ,'get_replaced_terms');

# test esearch data with History
$eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'esearch',
    -file       => test_input_file('eutils','esearch2.xml'));

is ($eutil->get_count, 534, 'get_count');
$history = $eutil->next_History;
isa_ok($history, 'Bio::Tools::EUtilities::HistoryI');
is($history->get_webenv,
   '00m7eJh8lyG3wiC2SE2hd7Im_w5o3z3q4_JK9-Rn266ix_eRXkjNOYQxHp@03F17619941CFD71_0005SID',
   'get_webenv');
is($history->get_query_key, 1, 'get_query_key');
is(join(',',$history->history),
   '00m7eJh8lyG3wiC2SE2hd7Im_w5o3z3q4_JK9-Rn266ix_eRXkjNOYQxHp@03F17619941CFD71_0005SID,1', 'history');

@ids2 = $eutil->get_ids;
is_deeply(\@ids2, \@ids, 'get_ids');
is($eutil->get_retstart, 0,'get_retstart');
is($eutil->get_retmax, 20,'get_retmax');
is($eutil->get_translation_from, 'Mus musculus','get_translation_from');
is($eutil->get_translation_to, '("Mus musculus"[Organism:__txid10090] OR Mus musculus)','get_translation_to');

# the database isn't carried into the parsers
is($eutil->get_db, undef, 'get_db');
is($eutil->get_database, undef, 'get_database');

# the term isn't carried into the parsers
is($eutil->get_term, undef,'get_term');

# espell only (should be undef)
is($eutil->get_corrected_query, undef ,'get_corrected_query');
is($eutil->get_replaced_terms, undef ,'get_replaced_terms');

my @qs = $eutil->get_GlobalQueries;
is(scalar(@qs), 0, 'get_GlobalQueries')
