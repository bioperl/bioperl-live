# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 13,
			   -requires_modules => [qw(IO::String
										LWP::UserAgent
										Bio::WebAgent
										HTML::HeadParser
										HTTP::Request::Common)],
			   -requires_networking => 1);
	
	use_ok('Bio::Tools::Analysis::DNA::ESEfinder');
	use_ok('Data::Dumper');
	use_ok('Bio::PrimarySeq');
}

#######all these tests work with 1ary seq########
my $seq = Bio::PrimarySeq->new(-id=>'bioperl',
                               -seq=>'atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt');
ok my $tool = Bio::Tools::Analysis::DNA::ESEfinder->new(-seq => $seq);

SKIP: {
	eval {$tool->run;};
	skip "Could not connect to ESEfinder server, skipping those tests", 9 if $@;
    ok my @res = $tool->result('Bio::SeqFeatureI');
	ok @res > 0;
    ok my $raw = $tool->result('');
    ok my $parsed = $tool->result('parsed');
    ok my $meta = $tool->result('all');
    is $parsed->[0][1], 41;
	
	SKIP: {
		test_skip(-tests => 3, -requires_module => 'Bio::Seq::Meta::Array');
		is $meta->{'seq'}, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt";
		is $meta->named_submeta_text('ESEfinder_SRp55', 1,2), "-3.221149 -1.602223";
		is $meta->seq, "atcgatgctatgcatgctatgggtgtgattcgatgcgactgttcatcgtagccccccccccccccctttt";
	}
}
