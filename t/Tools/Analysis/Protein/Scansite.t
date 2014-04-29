# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests            => 14,
               -requires_modules => [qw(IO::String
                                        LWP::UserAgent
                                        Data::Stag)]);

    use_ok('Bio::Tools::Analysis::Protein::Scansite');
    use_ok('Bio::SeqIO');
    use_ok('Bio::WebAgent');
}

my $verbose = test_debug();

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);


my $seqio=Bio::SeqIO->new( -verbose => $verbose,
                  -format => 'swiss',
                  -file   => test_input_file('swiss.dat'));

my $seq = $seqio->next_seq();
ok $tool = Bio::Tools::Analysis::Protein::Scansite->new( 
					-seq=>$seq->primary_seq);
ok $tool->stringency('Low');
is $tool->stringency(), 'Low';
is $tool->protein_id(), $tool->seq->display_id();

SKIP: {
	test_skip(-tests => 6, -requires_networking => 2);
	ok $tool->run();
	skip "Something wrong with server? Terminated by error, skipping tests", 5 if $tool->status eq 'TERMINATED_BY_ERROR';
	ok my $raw = $tool->result('');
	print $raw if $verbose;
	ok my $parsed = $tool->result('parsed');
	is $parsed->[0]{'site'}, 'T101';
	ok my @res = $tool->result('Bio::SeqFeatureI');
	is $res[0]->start, 101;
}
