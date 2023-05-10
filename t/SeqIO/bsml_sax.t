# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use Bio::Root::Test;
	
	test_begin(-tests => 15,
			   -requires_modules => [qw(XML::SAX
									    XML::SAX::Writer
									    XML::SAX::PurePerl
										XML::SAX::Base)]);
    
	use_ok('Bio::SeqIO');
}

my $verbose = test_debug();

# Using lax parser XML::SAX::PurePerl due to external DTD error (404 for URL,
# invalid redeclaration of predefined entity) instead of others such as
# XML::LibXML::SAX. See GH#376.
local $XML::SAX::ParserPackage = 'XML::SAX::PurePerl';
my $str = Bio::SeqIO->new(-format => 'bsml_sax',
			  -verbose => $verbose,
			  -file => test_input_file('U83300.bsml'));
my $seq = $str->next_seq;
isa_ok($seq, 'Bio::Seq::RichSeqI');
my @refs = $seq->annotation->get_Annotations('reference');
is(@refs, 2);
is($seq->display_id,'MIVN83300');
is($seq->molecule ,'dna');
ok(! $seq->is_circular);
is($seq->get_dates,2);
is($seq->accession_number, 'U83300');
is($seq->seq_version,1);
my @feats = $seq->get_SeqFeatures;
is(@feats, 2);
is($feats[1]->start, 1);
is($feats[1]->end, 946);
is($feats[1]->get_tag_values('db_xref'), 3);
is($seq->annotation->get_Annotations('reference'),2);
is($seq->annotation->get_Annotations('dblink'),2);
