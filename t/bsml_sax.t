# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($error $NUMTESTS);
BEGIN {
	$NUMTESTS = 16;
	$error = 0;
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if ( $@ ) {
		use lib 't\lib';
    }
    use Test::More;
	
	eval {
		require XML::SAX;
		require XML::SAX::Writer;
		require XML::SAX::Base;
	};
	if ($@) {
		plan skip_all => 'XML::SAX::Base or XML::SAX or XML::SAX::Writer not found - skipping bsml_sax tests';
	}
	else {
		plan tests => $NUMTESTS;
	}
}

use_ok('Bio::SeqIO');
use_ok('Bio::Root::IO');

my $verbose = $ENV{'BIOPERLDEBUG'};

my $str = Bio::SeqIO->new(-format => 'bsml_sax',
			  -verbose => $verbose,
			  -file => Bio::Root::IO->catfile
			  (qw(t data U83300.bsml) ));
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
