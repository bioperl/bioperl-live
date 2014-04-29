# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests               => 15,
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
					-seq=>$seq->primary_seq);
ok $tool->compartment(['golgi', 'er']);
ok my $cmp = $tool->compartment();
is $cmp->[1], 'GO:0005783';
ok $tool->species(9606);
is $tool->species, 9606;

ok $tool->run ();
exit if $tool->status eq 'TERMINATED_BY_ERROR';
ok my $raw = $tool->result('');
print $raw if $verbose;
ok my $parsed = $tool->result('parsed');
is $parsed->{'CLV_NDR_NDR_1'}{'locus'}[0], '54-56';
ok my @res = $tool->result('Bio::SeqFeatureI');
