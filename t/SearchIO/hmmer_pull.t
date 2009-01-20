# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_hmmer_pull.t 14984 2008-11-11 18:39:20Z sendu $

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 290);
	
	use_ok('Bio::SearchIO');
}

my $searchio = Bio::SearchIO->new(-format => 'hmmer_pull', -file => test_input_file('hmmpfam_fake.out'), -verbose => -1);
my @data = ([qw(roa1_drome roa2_drome)], [2, 1], [1, 2], [2, 1]);
while (my $result = $searchio->next_result) {
    is ref($result), 'Bio::Search::Result::HmmpfamResult';
    is $result->algorithm, 'HMMPFAM';
    is $result->algorithm_version, '2.1.1';
    is $result->hmm_name, 'pfam';
    is $result->hmm_file, $result->hmm_name;
    is $result->database_name, $result->hmm_name;
    is $result->sequence_file, '/home/birney/src/wise2/example/road.pep';
    is $result->sequence_database, $result->sequence_file;
    is $result->query_name, shift @{$data[0]};
    is $result->num_hits(), shift @{$data[1]};
    is $result->no_hits_found, 0;
    
	is $result->query_accession, '';
    is $result->query_description, '';
    ok ! $result->query_length;
    ok ! $result->database_letters;
    ok ! $result->database_entries;
    is $result->algorithm_reference, '';
    is $result->get_parameter('test'), undef;
    is $result->available_parameters, undef;
    is $result->get_statistic('test'), undef;
    is $result->available_statistics, undef;
    
    my @orig_order = $result->hits;
    is @orig_order, shift @{$data[3]};
	if (@orig_order > 1) {
		isnt $orig_order[0]->name, $orig_order[1]->name;
		$result->sort_hits(sub{$Bio::Search::Result::HmmpfamResult::a->[2]
													<=> 
							   $Bio::Search::Result::HmmpfamResult::b->[2]});
		my @hits = $result->hits;
		is @hits, @orig_order;
		is $hits[0]->name, $orig_order[1]->name;
		$result->sort_hits(sub{$Bio::Search::Result::HmmpfamResult::b->[4]
													<=> 
							   $Bio::Search::Result::HmmpfamResult::a->[4]});
	}
    
    my @hit_data = ([qw(SEED TEST)], [146.1, "5.0"], [6.3e-40, 7.2], [2, 1], [77, undef], [2, 0], [1, 2],
                    ["33 34 36 38 43 45 47 48 51 53 55 57 58 65 68 71 73 74 76 88 98 99 124 125 126 127 129 132 135 140 142 145 146 148 149 151 153 154 156 157 158 159 160 161 164 165 166 167 168 169 170 178 187 189 194", ''],
                    ["1 2 3 4 6 9 11 12 13 15 16 17 19 21 22 23 25 26 28 30 31 33 39 40 41 42 43 44 46 47 48 49 50 51 52 60 61 70 72 73 77", ''],
                    ["1-6 8-13 15-23 25-33 39-56 58-63 67-77", '']);
    while (defined(my $hit = $result->next_model)) {
        is ref($hit), 'Bio::Search::Hit::HmmpfamHit';
        is $hit->name, shift @{$hit_data[0]};
        is $hit->raw_score, shift @{$hit_data[1]};
        is $hit->score, $hit->raw_score;
        float_is $hit->significance, shift @{$hit_data[2]};
        float_is $hit->p, $hit->significance;
        is $hit->num_hsps, shift @{$hit_data[3]};
        is $hit->n, $hit->num_hsps;
        is $hit->algorithm, $result->algorithm;
        is $hit->overlap, 0;
        is $hit->rank, shift @{$hit_data[6]};
        is $hit->tiled_hsps, 0;
        is $hit->strand('query'), 1;
        is $hit->strand('hit'), 1;
        my @strands = $hit->strand;
        is "@strands", "1 1";
        
        is $hit->description, undef;
        is $hit->accession, undef;
        ok ! $hit->locus;
        ok ! $hit->bits;
        ok ! $result->logical_length('query');
        ok ! $result->frame;
        is $hit->each_accession_number, undef;
        
        is $hit->length, shift @{$hit_data[4]};
        is $hit->logical_length('hit'), $hit->length;
		
		if ($result->query_name eq 'roa1_drome') {
			my @inds = $hit->seq_inds('query', 'identical');
			is "@inds", shift @{$hit_data[7]};
			@inds = $hit->seq_inds('hit', 'identical');
			is "@inds", shift @{$hit_data[8]};
			@inds = $hit->seq_inds('hit', 'conserved', 1);
			is "@inds", shift @{$hit_data[9]};
		}
		
		if ($hit->name eq 'SEED') {
			my $best = $hit->hsp('best');
			float_is($best->evalue, 1.1e-18);
			my $worst = $hit->hsp('worst');
			float_is($worst->evalue, 2.2e-17);
			is $hit->start('query'), 33;
			is $hit->start('hit'), 1;
			is $hit->end('query'), 194;
			is $hit->end('hit'), 77;
			my @range = $hit->range('query');
			is "@range", '33 194';
			@range = $hit->range('hit');
			is "@range", '1 77';
			
			if ($hit->query_name eq 'roa1_drome') {
				is $hit->length_aln('query'),142;
				is $hit->length_aln('hit'), 77;
				is $hit->gaps('total'), 14;
				is $hit->gaps('query'), 13;
				is $hit->gaps('hit'), 1;
				is $hit->matches('id'), 41;
				is $hit->matches('cons'), 24;
				is $hit->frac_identical, 0.387;
				is $hit->frac_conserved, 0.169;
				ok ! $hit->frac_aligned_query;
				is $hit->frac_aligned_hit, '1.00';
				is $hit->num_unaligned_hit, 1;
				is $hit->num_unaligned_query, 13;
			}
		}
        
        my @hsps = $hit->hsps;
        is @hsps, shift @{$hit_data[5]};
        
        my @hsp_data = ([1, 1], [77, 77], [33, 124], [103, 194], [71.2, 75.5], [2.2e-17, 1.1e-18],
                        ['LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP',
                         'LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV'],
						[7, 6],
						['lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG.kelggrklrv',
                         'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv'],
                        ['lf+g+L + +t+e Lk++F+k G iv++ +++D     + t++s+Gf+F+++  ++  + A +    +++++gr+++ ',
                         'lfVg L  d +e+ ++d+F++fG iv+i+iv+D     ketgk +GfaFVeF++++ ++k +     ++l+g+ + v'],
						[1, 0], [8, 6], [1, 2], ['33 103', '124 194'], [78, 77], [22, 33], [33, 23],
						['0.3099', '0.4648'], ['0.2857', '0.4286'], ['0.2821', '0.4286']);
        
        while (defined(my $hsp = $hit->next_domain)) {
            is ref($hsp), 'Bio::Search::HSP::HmmpfamHSP';
            is $hsp->hit->start, shift @{$hsp_data[0]};
            is $hsp->hit->end, shift @{$hsp_data[1]};
            is $hsp->query->start, shift @{$hsp_data[2]};
            is $hsp->query->end, shift @{$hsp_data[3]};
			is $hsp->start('hit'), $hsp->hit->start;
			is $hsp->end('hit'),$hsp->hit->end;
			is $hsp->start('query'), $hsp->query->start;
			is $hsp->end('query'), $hsp->query->end;
			is $hsp->strand('hit'), 1;
			is $hsp->strand('query'), 1;
            is $hsp->score, shift @{$hsp_data[4]};
			ok ! $hsp->bits;
            float_is($hsp->evalue, shift @{$hsp_data[5]});
			ok ! $hsp->pvalue;
			float_is($hsp->significance, $hsp->evalue);
			is $hsp->algorithm, $result->algorithm;
			is $hsp->rank, shift @{$hsp_data[12]};
			my @range = $hsp->range;
			is "@range", shift @{$hsp_data[13]};
			is $hsp->n, $hit->num_hsps;
			is $hsp->length('query'), 71;
			is $hsp->length('hit'), 77;
			my $locseq = $hsp->seq('hit');
			
			if ($result->query_name eq 'roa1_drome') {
				is ref($locseq), 'Bio::LocatableSeq';
				my $aln = $hsp->get_aln('hit');
				is ref($aln), 'Bio::SimpleAlign';
				is $hsp->query_string, shift @{$hsp_data[6]};
				is $hsp->gaps('query'), shift @{$hsp_data[7]};
				is $hsp->gaps('hit'), shift @{$hsp_data[10]};
				is $hsp->gaps('total'), shift @{$hsp_data[11]};
				is $hsp->hit_string, shift @{$hsp_data[8]};
				is $hsp->homology_string, shift @{$hsp_data[9]};
				is $hsp->seq_str('hit'), $hsp->hit_string;
				is $hsp->seq_str('query'), $hsp->query_string;
				is $hsp->seq_str('homology'), $hsp->homology_string;
				is length($hsp->homology_string), length($hsp->hit_string);
				is length($hsp->query_string), length($hsp->homology_string);
				is $hsp->length('total'), shift @{$hsp_data[14]};
				is $hsp->hsp_length, $hsp->length('total');
				is $hsp->num_identical, shift @{$hsp_data[15]};
				is $hsp->num_conserved, shift @{$hsp_data[16]};
				is $hsp->frac_identical('query'), shift @{$hsp_data[17]};
				is $hsp->frac_identical('hit'), shift @{$hsp_data[18]};
				is $hsp->frac_identical('total'), shift @{$hsp_data[19]};
			}
        }
    }
}

is $searchio->result_count, 2;

# bug revealed by bug 2632 - CS lines were already ignored, but we couldn't
# parse alignments when HSPs weren't in simple order!!
$searchio = Bio::SearchIO->new(-format => 'hmmer_pull', -file => test_input_file('hmmpfam_cs.out'), -verbose => 1);
my $result = $searchio->next_result;
my $hit = $result->next_hit;
my $hsp = $hit->next_hsp;
is $hsp->seq_str, "IPPLLAVGAVHHHLINKGLRQEASILV";

# and another bug revealed: we don't always know the hit length, and
# shouldn't complain about that with a warning
is $hsp->hit->seqlength, 412;

my $count = 0;
while (my $hit = $result->next_hit) {
    $count++;
    next if $count < 6;
    last if $count > 6;
	my $hsp = $hit->next_hsp;
    ok ! $hsp->hit->seqlength;
    #*** not sure how to test for the lack of a warning though...
	# Maybe run an eval with verbose set to 2, then make sure $@ is undef? --cjfields
}
