# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 10,
			   -requires_modules => [qw(IO::String LWP::UserAgent)],
			   -requires_networking => 1);
	
	use_ok 'Bio::Tools::Analysis::Protein::Mitoprot';
	use_ok 'Bio::PrimarySeq';
	use_ok 'Bio::WebAgent';
}

my $verbose = test_debug();

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

my $seq = Bio::PrimarySeq->new(-seq => 'MSADQRWRQDSQDSFGDSFDGDSFFGSDFDGDS'.
                               'DFGSDFGSDGDFGSDFGDSFGDGFSDRSRQDQRS',
                               -display_id => 'test2');

ok $tool = Bio::Tools::Analysis::Protein::Mitoprot->new( -seq=>$seq);
SKIP: {
	ok $tool->run();
	skip('Server terminated with an error, skipping tests', 4) if $tool->status eq 'TERMINATED_BY_ERROR';
	ok my $raw = $tool->result('');
	ok my $parsed = $tool->result('parsed');
	is ($parsed->{'charge'}, -13);
	ok my @res = $tool->result('Bio::SeqFeatureI');
}
