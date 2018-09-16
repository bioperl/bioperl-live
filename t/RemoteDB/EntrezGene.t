# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 8,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent
										HTTP::Request::Common
                                        Bio::ASN1::EntrezGene)
                                        ],
			   -requires_networking => 1);
	use_ok('Bio::DB::EntrezGene');
}

my ($gb, $seq, $seqio);

#
# Bio::DB::EntrezGene
#
SKIP: {
	test_skip(-tests => 7);
    ok $gb = Bio::DB::EntrezGene->new(-retrievaltype => 'tempfile', -delay => 0);
    eval {$seqio = $gb->get_Stream_by_id([2,3064]);};
    skip "Couldn't connect to Entrez with Bio::DB::EntrezGene. Skipping those tests", 6 if $@;
    $seq = $seqio->next_seq;
    is $seq->display_id, "A2M";
    is $seq->accession_number, 2;
    $seq = $seqio->next_seq;
    is $seq->display_id, "HTT";
    is $seq->accession_number, 3064;
    eval {$seq = $gb->get_Seq_by_id(6099);};
    skip "Couldn't connect to Entrez with Bio::DB::EntrezGene. Skipping those tests", 2 if $@;
    is $seq->display_id, "RP";
    is $seq->accession_number, 6099;
}
