# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan test => 287;
}

use Bio::Root::IO;
use Bio::SearchIO;

ok(1);

my $searchio = new Bio::SearchIO(-format => 'hmmer_pull', -file => Bio::Root::IO->catfile("t","data","hmmpfam_fake.out"), -verbose => -1);
my @data = ([qw(roa1_drome roa2_drome)], [2, 1], [1, 2], [2, 1]);
while (my $result = $searchio->next_result) {
    ok ref($result), 'Bio::Search::Result::HmmpfamResult';
    ok $result->algorithm, 'HMMPFAM';
    ok $result->algorithm_version, '2.1.1';
    ok $result->hmm_name, 'pfam';
    ok $result->hmm_file, $result->hmm_name;
    ok $result->database_name, $result->hmm_name;
    ok $result->sequence_file, '/home/birney/src/wise2/example/road.pep';
    ok $result->sequence_database, $result->sequence_file;
    ok $result->query_name, shift @{$data[0]};
    ok $result->num_hits(), shift @{$data[1]};
    ok $result->no_hits_found, 0;
    
	ok $result->query_accession, '';
    ok $result->query_description, '';
    ok ! $result->query_length;
    ok ! $result->database_letters;
    ok ! $result->database_entries;
    ok $result->algorithm_reference, '';
    ok $result->get_parameter('test'), undef;
    ok $result->available_parameters, undef;
    ok $result->get_statistic('test'), undef;
    ok $result->available_statistics, undef;
    
    my @orig_order = $result->hits;
    ok @orig_order, shift @{$data[3]};
	if (@orig_order > 1) {
		ok $orig_order[0]->name ne $orig_order[1]->name;
		$result->sort_hits(sub{$Bio::Search::Result::HmmpfamResult::a->[2]
													<=> 
							   $Bio::Search::Result::HmmpfamResult::b->[2]});
		my @hits = $result->hits;
		ok @hits, @orig_order;
		ok $hits[0]->name, $orig_order[1]->name;
		$result->sort_hits(sub{$Bio::Search::Result::HmmpfamResult::b->[4]
													<=> 
							   $Bio::Search::Result::HmmpfamResult::a->[4]});
	}
    
    my @hit_data = ([qw(SEED TEST)], [146.1, "5.0"], [6.3e-40, 7.2], [2, 1], [77, undef], [2, 0], [1, 2],
                    ["33 34 36 38 43 45 47 48 51 53 55 57 58 65 68 71 73 74 76 88 98 99 124 125 126 127 129 132 135 140 142 145 146 148 149 151 153 154 156 157 158 159 160 161 164 165 166 167 168 169 170 178 187 189 194", ''],
                    ["1 2 3 4 6 9 11 12 13 15 16 17 19 21 22 23 25 26 28 30 31 33 39 40 41 42 43 44 46 47 48 49 50 51 52 60 61 70 72 73 77", ''],
                    ["1-6 8-13 15-23 25-33 39-56 58-63 67-77", '']);
    while (defined(my $hit = $result->next_model)) {
        ok ref($hit), 'Bio::Search::Hit::HmmpfamHit';
        ok $hit->name, shift @{$hit_data[0]};
        ok $hit->raw_score, shift @{$hit_data[1]};
        ok $hit->score, $hit->raw_score;
        ok $hit->significance, shift @{$hit_data[2]};
        ok $hit->p, $hit->significance;
        ok $hit->num_hsps, shift @{$hit_data[3]};
        ok $hit->n, $hit->num_hsps;
        ok $hit->algorithm, $result->algorithm;
        ok $hit->overlap, 0;
        ok $hit->rank, shift @{$hit_data[6]};
        ok $hit->tiled_hsps, 0;
        ok $hit->strand('query'), 1;
        ok $hit->strand('hit'), 1;
        my @strands = $hit->strand;
        ok "@strands", "1 1";
        
        ok $hit->description, undef;
        ok $hit->accession, undef;
        ok ! $hit->locus;
        ok ! $hit->bits;
        ok ! $result->logical_length('query');
        ok ! $result->frame;
        ok $hit->each_accession_number, undef;
        
        ok $hit->length, shift @{$hit_data[4]};
        ok $hit->logical_length('hit'), $hit->length;
		
		if ($result->query_name eq 'roa1_drome') {
			my @inds = $hit->seq_inds('query', 'identical');
			ok "@inds", shift @{$hit_data[7]};
			@inds = $hit->seq_inds('hit', 'identical');
			ok "@inds", shift @{$hit_data[8]};
			@inds = $hit->seq_inds('hit', 'conserved', 1);
			ok "@inds", shift @{$hit_data[9]};
		}
		
		if ($hit->name eq 'SEED') {
			my $best = $hit->hsp('best');
			ok $best->evalue, 1.1e-18;
			my $worst = $hit->hsp('worst');
			ok $worst->evalue, 2.2e-17;
			ok $hit->start('query'), 33;
			ok $hit->start('hit'), 1;
			ok $hit->end('query'), 194;
			ok $hit->end('hit'), 77;
			my @range = $hit->range('query');
			ok "@range", '33 194';
			@range = $hit->range('hit');
			ok "@range", '1 77';
			
			if ($hit->query_name eq 'roa1_drome') {
				ok $hit->length_aln('query'),142;
				ok $hit->length_aln('hit'), 77;
				ok $hit->gaps('total'), 14;
				ok $hit->gaps('query'), 13;
				ok $hit->gaps('hit'), 1;
				ok $hit->matches('id'), 41;
				ok $hit->matches('cons'), 24;
				ok $hit->frac_identical, 0.387;
				ok $hit->frac_conserved, 0.169;
				ok ! $hit->frac_aligned_query;
				ok $hit->frac_aligned_hit, '1.00';
				ok $hit->num_unaligned_hit, 1;
				ok $hit->num_unaligned_query, 13;
			}
		}
        
        my @hsps = $hit->hsps;
        ok @hsps, shift @{$hit_data[5]};
        
        my @hsp_data = ([1, 1], [77, 77], [33, 124], [103, 194], [71.2, 75.5], [2.2e-17, 1.1e-18],
                        ['LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP',
                         'LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV'],
						[7, 6],
						['lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG.kelggrklrv',
                         'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv'],
                        ['lf+g+L + +t+e Lk++F+k G iv++ +++D     + t++s+Gf+F+++  ++  + A +    +++++gr+++ ',
                         'lfVg L  d +e+ ++d+F++fG iv+i+iv+D     ketgk +GfaFVeF++++ ++k +     ++l+g+ + v'],
						[1, 0], [8, 6], [1, 2], ['33 103', '124 194'], [78, 77], [22, 33], [33, 23],
						['0.310', '0.465'], ['0.286', '0.429'], ['0.282', '0.429']);
        
        while (defined(my $hsp = $hit->next_domain)) {
            ok ref($hsp), 'Bio::Search::HSP::HmmpfamHSP';
            ok $hsp->hit->start, shift @{$hsp_data[0]};
            ok $hsp->hit->end, shift @{$hsp_data[1]};
            ok $hsp->query->start, shift @{$hsp_data[2]};
            ok $hsp->query->end, shift @{$hsp_data[3]};
			ok $hsp->start('hit'), $hsp->hit->start;
			ok $hsp->end('hit'),$hsp->hit->end;
			ok $hsp->start('query'), $hsp->query->start;
			ok $hsp->end('query'), $hsp->query->end;
			ok $hsp->strand('hit'), 1;
			ok $hsp->strand('query'), 1;
            ok $hsp->score, shift @{$hsp_data[4]};
			ok ! $hsp->bits;
            ok $hsp->evalue, shift @{$hsp_data[5]};
			ok ! $hsp->pvalue;
			ok $hsp->significance, $hsp->evalue;
			ok $hsp->algorithm, $result->algorithm;
			ok $hsp->rank, shift @{$hsp_data[12]};
			my @range = $hsp->range;
			ok "@range", shift @{$hsp_data[13]};
			ok $hsp->n, $hit->num_hsps;
			ok $hsp->length('query'), 71;
			ok $hsp->length('hit'), 77;
			my $locseq = $hsp->seq('hit');
			
			if ($result->query_name eq 'roa1_drome') {
				ok ref($locseq), 'Bio::LocatableSeq';
				my $aln = $hsp->get_aln('hit');
				ok ref($aln), 'Bio::SimpleAlign';
				ok $hsp->query_string, shift @{$hsp_data[6]};
				ok $hsp->gaps('query'), shift @{$hsp_data[7]};
				ok $hsp->gaps('hit'), shift @{$hsp_data[10]};
				ok $hsp->gaps('total'), shift @{$hsp_data[11]};
				ok $hsp->hit_string, shift @{$hsp_data[8]};
				ok $hsp->homology_string, shift @{$hsp_data[9]};
				ok $hsp->seq_str('hit'), $hsp->hit_string;
				ok $hsp->seq_str('query'), $hsp->query_string;
				ok $hsp->seq_str('homology'), $hsp->homology_string;
				ok length($hsp->homology_string), length($hsp->hit_string);
				ok length($hsp->query_string), length($hsp->homology_string);
				ok $hsp->length('total'), shift @{$hsp_data[14]};
				ok $hsp->hsp_length, $hsp->length('total');
				ok $hsp->num_identical, shift @{$hsp_data[15]};
				ok $hsp->num_conserved, shift @{$hsp_data[16]};
				ok $hsp->frac_identical('query'), shift @{$hsp_data[17]};
				ok $hsp->frac_identical('hit'), shift @{$hsp_data[18]};
				ok $hsp->frac_identical('total'), shift @{$hsp_data[19]};
			}
        }
    }
}

ok $searchio->result_count, 2;
