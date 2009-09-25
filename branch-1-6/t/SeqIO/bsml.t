# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 16,
			   -requires_modules => ['XML::DOM']
			  );
	use_ok('XML::DOM');
	use_ok('Bio::SeqIO::bsml');
}

my $verbose = test_debug();

my $str = Bio::SeqIO->new(-format => 'bsml',
			  -verbose => $verbose,
			  -file => test_input_file('U83300.bsml'));
my $seq = $str->next_seq;
isa_ok($seq, 'Bio::Seq::RichSeqI');
my @refs = $seq->annotation->get_Annotations('reference');
is(@refs, 2, 'got correct number of refs');
is($seq->display_id, 'MIVN83300', 'display_id');
is($seq->molecule, 'DNA', 'molecule');
ok(! $seq->is_circular, 'is_circular');
is($seq->get_dates, 2, 'dates');
is($seq->accession_number, 'U83300', 'accession_number');
is($seq->seq_version, 1, 'seq_version');
my @feats = $seq->get_SeqFeatures;
is(@feats, 2, 'got correct number of SeqFeatures');
is($feats[1]->start, 1, 'feature start');
is($feats[1]->end, 946, 'feature end');
is($feats[1]->get_tag_values('db_xref'), 3, 'get_tag_values db_xref');
is($seq->annotation->get_Annotations('reference'), 2, 'get_Annotations reference');
is($seq->annotation->get_Annotations('dblink'), 2, 'get_Annotations dblink');
