# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 16,
               -requires_modules => [qw(IO::String LWP::UserAgent)]);
	
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Tools::Analysis::Protein::Sopma');
}

my $verbose = test_debug();

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

my $seq = Bio::PrimarySeq->new(
  -seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
  -display_id => 'test2');
ok $tool = Bio::Tools::Analysis::Protein::Sopma->new( -seq=>$seq,
													  #-verbose => $verbose,
                                                      -window_width => 15);

SKIP: {
	test_skip(-tests => 12, -requires_networking => 1);
	
	ok $tool->run();
	skip "Tool was terminated by some error: problem connecting to server?", 11 if $tool->status eq 'TERMINATED_BY_ERROR';
	
	ok my $raw = $tool->result('');
	ok my $parsed = $tool->result('parsed');
	is ($parsed->[0]{'helix'}, '102');
	ok my @res = $tool->result('Bio::SeqFeatureI');
	ok my $meta = $tool->result('meta', "ww15");

	ok $tool->window_width(21);
	ok $tool->clear();
	ok $tool->run;
	ok my $meta2 = $tool->result('meta', "ww21");
	
	SKIP: {
		test_skip(-tests => 2, -requires_module => 'Bio::Seq::Meta::Array');
		is $meta->named_submeta_text('Sopma_helix|ww15',1,2), '102 195';
		is $meta->seq, 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS';
	}
}
