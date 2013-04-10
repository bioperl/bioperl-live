# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 26,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent)],
			   -requires_networking => 1);
	
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Tools::Analysis::Protein::Domcut');
}

my $verbose = test_debug();

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

SKIP: {
	######## test using PrimarySeq object ##############
	my $seq = Bio::PrimarySeq->new(-seq        => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQPPPPPPPPPPPPPDQRS',
								   -display_id => 'test2');
	
	ok $tool = Bio::Tools::Analysis::Protein::Domcut->new( -seq=>$seq);
	ok $tool->run ();
	if ($tool->status eq 'TERMINATED_BY_ERROR') {
		skip('Problem with DomCut run, check status', 21);
	}
	
	ok my $raw    = $tool->result('');
	ok my $parsed = $tool->result('parsed');
	is ($parsed->[23]{'score'}, '-0.209');
	my @res       = $tool->result('Bio::SeqFeatureI');
	if (scalar @res > 0) {
		ok 1;
	} else {
		skip('No network access - could not connect to Domcut server', 18);
	}
	ok my $meta = $tool->result('meta');
	
	SKIP: {
		test_skip(-tests => 2, -requires_module => 'Bio::Seq::Meta::Array');
		is($meta->named_submeta_text('Domcut', 1,2), "0.068 0.053");
		is ($meta->seq, "MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQPPPPPPPPPPPPPDQRS");
	}
	
	########## test using Bio::Seq object ##############
	ok my $tool2 = Bio::WebAgent->new(-verbose =>$verbose);
	
	ok my $seq2  = Bio::Seq->new(-seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
				 -display_id => 'test2');
	
	ok $tool2 = Bio::Tools::Analysis::Protein::Domcut->new( -seq=>$seq2->primary_seq);
	ok $tool2->run ();
	
	@res = $tool2->result('Bio::SeqFeatureI');
	
	if (scalar @res > 0) {
		ok 1;
	} else {
		skip('No network access - could not connect to Domcut server', 10);
	}
	
	ok my $parsed2 = $tool2->result('parsed');
	is ($parsed2->[23]{'score'}, '-0.209');
	
	ok my $meta2 = $tool2->result('meta');
	
	is($meta2->named_submeta_text('Domcut', 1,2), "0.068 0.053");
	is ($meta2->seq, "MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS");
	
	ok my $seq4 = Bio::Seq->new();
	ok $seq2->primary_seq($meta2);
	for (@res) {
		ok $seq2->add_SeqFeature($_);
	}
	ok $seq2->primary_seq->named_submeta_text('Domcut', 1,2);
}
