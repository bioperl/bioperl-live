# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_hmmer3.t 14989 2008-11-11 19:52:02Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 10);

	use_ok('Bio::SearchIO');
}

my $searchio = Bio::SearchIO->new(-format  => 'hmmer3',
                                  -file    => test_input_file('hmmscan.out'),
                                  -verbose => 1);

while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::hmmer3Result', 'Check for the correct reference type');
    is($result->algorithm, 'HMMSCAN', 'Check algorithm');
    is($result->algorithm_version, '3.0', 'Check algorithm version');
    is($result->hmm_name, '/data/biodata/HMMerDB/Pfam.hmm', 'Check hmm_name');
    is($result->sequence_file, 'BA000019.orf1.fasta', 'Check sequence_file');
    is($result->query_name, 'BA000019.orf1', 'Check query_name');
    is($result->query_length, '198', 'Check query_length');
    is($result->query_description, '', 'Check query_description');
    is($result->num_hits(), 1, 'Check num_hits');
}
