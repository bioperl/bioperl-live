# -*-Perl-*- Test Harness script for Bioperl
# $Id: epost.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 17,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'epost',
    -file       => test_input_file('eutils','epost.xml'));

# all parsers and data objects implement eutil() and datatype() (generally for
# debuggin purposes, but others may find them useful)
isa_ok($eutil, 'Bio::Tools::EUtilities::EUtilDataI');
isa_ok($eutil, 'Bio::Tools::EUtilities::Query');
is($eutil->eutil, 'epost', 'eutil');
is($eutil->datatype, 'query', 'datatype');
my $history = $eutil->next_History;
isa_ok($history, 'Bio::Tools::EUtilities::HistoryI');
isa_ok($history, 'Bio::Tools::EUtilities::EUtilDataI');
is($history->eutil, 'epost', 'eutil');
is($history->datatype, 'history', 'eutil');

# simple epost does not have anything other than the webenv/query_key
is($history->get_webenv,
   '0rACq8_iP87yHkqqm0SBaU38LzWLHIUd-J4QozMr31bh_XO5KAxLr5Q0o2e@03ED1E11941B69F1_0100SID',
   'get_webenv');
is($history->get_query_key, 1, 'get_query_key');
is(join(',',$history->history),
   '0rACq8_iP87yHkqqm0SBaU38LzWLHIUd-J4QozMr31bh_XO5KAxLr5Q0o2e@03ED1E11941B69F1_0100SID,1', 'history');
is($eutil->get_database, undef, 'get_database');
is($eutil->get_ids, undef, 'get_ids');

my @ids = qw(1621261 89318838 68536103 20807972 730439);

# add Parameters
my $pb = Bio::Tools::EUtilities::EUtilParameters->new(-eutil => 'epost',
									   -db => 'protein',
									   -id => \@ids);

$eutil->parameter_base($pb);

is($eutil->get_database, 'protein', 'get_database');
my @ids2 = $eutil->get_ids;
is_deeply(\@ids2, \@ids, 'get_ids');
