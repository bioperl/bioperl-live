# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_hmmer3.t 14989 2008-11-11 19:52:02Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 116);

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
        is($hit->description, 'Periplasmic binding protein', 'Check for hit description');
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

$searchio = Bio::SearchIO->new(-format  => 'hmmer3',
                               -file    => test_input_file('hmmsearch3.out'),
                               -verbose => 1);
while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::hmmer3Result', 'Check for the correct result reference type');
    is($result->algorithm, 'HMMSEARCH', 'Check algorithm');
    is($result->algorithm_version, '3.0', 'Check algorithm version');
    is($result->hmm_name, 'Kv9.hmm', 'Check hmm_name');
    is($result->sequence_file, '/home/pboutet/Desktop/databases/nr_May26', 'Check sequence_file');
    is($result->query_name, 'Kv9', 'Check query_name');
    is($result->query_length, '481', 'Check query_length');
    is($result->query_description, '', 'Check query_description');
    is($result->num_hits(), 2, 'Check num_hits');
    while( my $hit = $result->next_model ) {
    }
}

$searchio = Bio::SearchIO->new(-format  => 'hmmer3',
                               -file    => test_input_file('hmmscan_multi_domain.out'),
                               -verbose => 1);

my @multi_hits = (
      ['PPC', 'Bacterial pre-peptidase C-terminal domain', '111.0', 3.1e-32, 6,
        [[84, 192, 4, 59, 0.5, 0.16], [311, 397, 12, 58, -0.6, 0.36],
         [470, 550, 1, 69, 71.3, 1.3e-23], [567, 626, 15, 25, -3.2, 2],
         [974, 1072, 13, 36, -1.1, 0.5], [1087, 1169, 1, 69, 54.4, 2.4e-18]]],
      ['HemolysinCabind', 'Hemolysin-type calcium-binding repeat (2 cop', '47.9', 4.7e-13, 3,
        [[1213, 1225, 2, 13, 5.9, 0.0026], [1231, 1248, 1, 18, 10.8, 6.8e-5],
         [1240, 1257, 4, 18, 11.4, 4.3e-05]]]
    );

while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::hmmer3Result', 'Check for the correct result reference type');
    is($result->algorithm, 'HMMSCAN', 'Check algorithm');
    is($result->algorithm_version, '3.0', 'Check algorithm version');
    is($result->hmm_name, '/data/biodata/HMMerDB/Pfam-A.hmm', 'Check hmm_name');
    is($result->sequence_file, 'BA000019.orf37.fasta', 'Check sequence_file');
    is($result->query_name, 'BA000019.orf37', 'Check query_name');
    is($result->query_length, '1418', 'Check query_length');
    is($result->query_description, '', 'Check query_description');
    is($result->num_hits(), 2, 'Check num_hits');
    my ($hsp,$hit);
    while( $hit = $result->next_model ) {
        my @expected = @{shift @multi_hits};
        is(ref($hit), 'Bio::Search::Hit::hmmer3Hit', 'Check for the correct hit reference type');
        is($hit->name, shift @expected, 'Check hit name');
        is($hit->description, shift @expected, 'Check for hit description');
        is($hit->raw_score, shift @expected, 'Check hit raw_score');
        float_is($hit->significance, shift @expected, 'Check hit significance');
        is($hit->num_hsps, shift @expected, 'Check num_hsps');
        my @hsp_list = @{shift @expected};
        while( defined( $hsp = $hit->next_domain ) ) {
            my @hsp_exp = @{shift @hsp_list};
            is(ref($hsp), 'Bio::Search::HSP::hmmer3HSP', 'Check for correct hsp reference type');
            is($hsp->hit->start, shift @hsp_exp, 'Check for hit envfrom value');
            is($hsp->hit->end, shift @hsp_exp, 'Check for hit env to value');
            is($hsp->query->start, shift @hsp_exp, 'Check for query hmmfrom value');
            is($hsp->query->end, shift @hsp_exp, 'Check for query hmm to value');
            is($hsp->score, shift @hsp_exp, 'Check for hsp score');
            float_is($hsp->evalue, shift @hsp_exp, 'Check for hsp c-Evalue');
        }
    }
}
