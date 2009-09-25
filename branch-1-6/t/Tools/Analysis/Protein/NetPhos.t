# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14,
			   -requires_modules => [qw(IO::String LWP::UserAgent)],
			   -requires_networking => 1);
	
	use_ok('Bio::Tools::Analysis::Protein::NetPhos');
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::WebAgent');
}

my $verbose = test_debug();

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

ok $tool->sleep;
is $tool->delay(1), 1;
ok $tool->sleep;
ok $tool->timeout(120); # LWP::UserAgent method
is $tool->url('http://a.b.c/'), 'http://a.b.c/';


my $seq = Bio::PrimarySeq->new(-id=>'bioperl',
							   -seq=>'ABCDEFGHIJKLLKJFHSAKNDJFPSINCSJNDSKNSN');

ok $tool = Bio::Tools::Analysis::Protein::NetPhos->new(-verbose =>$verbose);
$tool->timeout(15);

ok $tool->run ( {seq=>$seq, threshold=>0.9} );
SKIP: {
	if ($tool->status eq 'TERMINATED_BY_ERROR') {
		skip "Running of the tool was terminated by an error, probably network/ NetPhos server error", 3;
	}
	my @res = $tool->result('Bio::SeqFeatureI');
	unless (@res) {
		skip "Didn't get any results from NetPhos server, probable network/server error", 3;
	}
	#new tests her in v 1.2
	ok my $raw = $tool->result('');
	ok my $parsed = $tool->result('parsed');
	is $parsed->[0][1], '0.934';
}
