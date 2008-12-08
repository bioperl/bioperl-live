# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 13,
			   -requires_modules => [qw(IO::String LWP::UserAgent)],
			   -requires_networking => 1);
	
	use_ok("Bio::Seq");
	use_ok("Bio::Tools::Analysis::Protein::GOR4");
}

my $seq = Bio::Seq->new(-seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
                        -display_id => 'test2');
ok my $tool = Bio::Tools::Analysis::Protein::GOR4->new(-seq=>$seq->primary_seq);

SKIP: {
    ok $tool->run();
	skip "Skipping tests since we got terminated by a server error", 9 if $tool->status eq 'TERMINATED_BY_ERROR';
    ok my $raw = $tool->result('');
    ok my $parsed = $tool->result('parsed');
    is $parsed->[0]{'coil'}, '999';
    my @res = sort {$a->start <=> $b->start} $tool->result('Bio::SeqFeatureI');
    if (scalar @res > 0) {
		ok 1;
    }
	else {
		skip 'No results - could not connect to GOR4 server?', 6;
    }
	is $res[0]->start, 1;
	is $res[0]->end, 43;
    ok my $meta = $tool->result('meta');
    
	test_skip(-tests => 2, -requires_module => 'Bio::Seq::Meta::Array');
	is $meta->named_submeta_text('GOR4_coil',1,2), '999 999';
	is $meta->seq, 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS';
}
