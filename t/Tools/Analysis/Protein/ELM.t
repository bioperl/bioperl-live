# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
BEGIN {
    use Bio::Root::Test;

    test_begin(-tests               => 16,
               -requires_modules    => [qw(IO::String
                                           LWP::UserAgent
                                           HTML::HeadParser
                                           Data::Stag)],
               -requires_networking => 1);

    use_ok('Bio::Tools::Analysis::Protein::ELM');
    use_ok('Bio::SeqIO');
    use_ok('Bio::WebAgent');
}

my $verbose = test_debug();

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

my $seqio=Bio::SeqIO->new( -verbose => $verbose,
                  -format => 'swiss',
                  -file   => test_input_file('swiss.dat'));

my $seq = $seqio->next_seq();
ok $tool = Bio::Tools::Analysis::Protein::ELM->new(
					-seq=>$seq->primary_seq), 'new object';
ok $tool->compartment(['golgi', 'er']), 'set compartment';
ok my $cmp = $tool->compartment(), 'get compartment';
is $cmp->[1], 'GO:0005783', 'check compartment';
ok $tool->species(9606), 'set species()';
is $tool->species, 9606, 'get species()';;

my $req_status = $tool->run();

ok $req_status, 'run';

my $status = $tool->status();

ok(defined($status), 'This returns something valid');

SKIP: {
    skip "Bad run() status, possible time out or error so skipping tests", 4 if !$req_status or $status eq 'TERMINATED_BY_ERROR';

    ok my $raw = $tool->result('');
    print $raw if $verbose;
    ok my $parsed = $tool->result('parsed');

    is $parsed->{'CLV_NRD_NRD_1'}{'locus'}[0], '54-56';
    ok my @res = $tool->result('Bio::SeqFeatureI');
};
