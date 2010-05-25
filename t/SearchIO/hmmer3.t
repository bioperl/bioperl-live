# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_hmmer3.t 14989 2008-11-11 19:52:02Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 22);

	use_ok('Bio::SearchIO');
}

my $searchio = Bio::SearchIO->new(-format  => 'hmmer3',
                                  -file    => test_input_file('hmmscan.out'),
                                  -verbose => 1);

while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::hmmer3Result', 'Check for the correct result reference type');
    is($result->algorithm, 'HMMSCAN', 'Check algorithm');
    is($result->algorithm_version, '3.0', 'Check algorithm version');
    is($result->hmm_name, '/data/biodata/HMMerDB/Pfam.hmm', 'Check hmm_name');
    is($result->sequence_file, 'BA000019.orf1.fasta', 'Check sequence_file');
    is($result->query_name, 'BA000019.orf1', 'Check query_name');
    is($result->query_length, '198', 'Check query_length');
    is($result->query_description, '', 'Check query_description');
    is($result->num_hits(), 1, 'Check num_hits');
    my ($hsp,$hit);
    if( $hit = $result->next_model ) {
        is(ref($hit), 'Bio::Search::Hit::hmmer3Hit', 'Check for the correct hit reference type');
        is($hit->name, 'Peripla_BP_2', 'Check hit name');
        is($hit->raw_score, '105.2', 'Check hit raw_score');
        float_is($hit->significance, 6e-30, 'Check hit significance');
        is($hit->num_hsps, 1, 'Check num_hsps');

        if( defined( $hsp = $hit->next_domain ) ) {
            is(ref($hsp), 'Bio::Search::HSP::hmmer3HSP', 'Check for correct hsp reference type');
            is($hsp->hit->start, 1, 'Check for hit envfrom value');
            is($hsp->hit->end, 175, 'Check for hit env to value');
            is($hsp->query->start, 59, 'Check for query hmmfrom value');
            is($hsp->query->end, 236, 'Check for query hmm to value');
            is($hsp->score, '105.0', 'Check for hsp score');
            float_is($hsp->evalue, 1.5e-33, 'Check for hsp c-Evalue');
        }
    }
}
