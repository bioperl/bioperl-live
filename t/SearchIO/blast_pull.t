# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_blast_pull.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 289);
	
	use_ok('Bio::SearchIO');
}


my $searchio = Bio::SearchIO->new(-format => 'blast_pull', -file => test_input_file('new_blastn.txt'));

my $result = $searchio->next_result;
is $result->database_name, 'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS,GSS,environmental samples or phase 0, 1 or 2 HTGS sequences)';
is $result->database_entries, 3742891;
is $result->database_letters, 16670205594;
is $result->algorithm, 'BLASTN';
is $result->algorithm_version, '2.2.13 [Nov-27-2005]';
is $result->query_name, 'pyrR,';
is $result->query_length, 558;
is $result->get_statistic('kappa'), 0.711;
is $result->get_statistic('kappa_gapped'), 0.711;
is $result->get_statistic('lambda'), 1.37;
is $result->get_statistic('lambda_gapped'), 1.37;
is $result->get_statistic('entropy'), 1.31;
is $result->get_statistic('entropy_gapped'), 1.31;
is $result->get_statistic('dbletters'), -509663586;
is $result->get_statistic('dbentries'), 3742891;
is $result->get_statistic('effectivespace'), 8935230198384;
is $result->get_parameter('matrix'), 'blastn matrix:1 -3';
is $result->get_parameter('gapopen'), 5;
is $result->get_parameter('gapext'), 2;
is $result->get_statistic('S2'), 60;
is $result->get_statistic('S2_bits'), 119.4;
float_is $result->get_parameter('expect'), 1e-23;
is $result->get_statistic('num_extensions'), 117843;

my @valid = ( [ 'gi|41400296|gb|AE016958.1|', 4829781, 'AE016958', '6e-059', 236],
	      [ 'gi|54013472|dbj|AP006618.1|', 6021225, 'AP006618', '4e-026', 127],
	      [ 'gi|57546753|dbj|BA000030.2|', 9025608, 'BA000030', '1e-023', 119]);
my $count = 0;
while (my $hit = $result->next_hit ) {
    my $d = shift @valid;

    is $hit->name, shift @$d;
    is $hit->length, shift @$d;
    is $hit->accession, shift @$d;
    float_is($hit->significance, shift @$d);
    is $hit->raw_score, shift @$d;

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is $hsp->query->start, 262;
            is $hsp->query->end, 552;
            is $hsp->hit->start, 1166897;
            is $hsp->hit->end, 1167187;
            is $hsp->length('hsp'), 291;
            is $hsp->start('hit'), $hsp->hit->start;
            is $hsp->end('query'), $hsp->query->end;
            is $hsp->strand('sbjct'), $hsp->subject->strand;# alias for hit
            float_is($hsp->evalue, 6e-59);
            is $hsp->score, 119;
            is $hsp->bits,236;	    	    
            is sprintf("%.1f",$hsp->percent_identity), 85.2;
            is sprintf("%.3f",$hsp->frac_identical('query')), 0.852;
            is sprintf("%.3f",$hsp->frac_identical('hit')), 0.852;
            is $hsp->gaps, 0;
            $hsps_left--;
        }
        is $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
is @valid, 0;

# descriptionless hit
$searchio = Bio::SearchIO->new('-format' => 'blast_pull',
							   '-file' => test_input_file('blast_no_hit_desc.txt'));

$result = $searchio->next_result;
my $hit = $result->next_hit;
is $hit->name, 'chr1';
is $hit->description, '';

# further (NCBI blastn/p) tests taken from SearchIO.t
$searchio = Bio::SearchIO->new('-format' => 'blast_pull',
							   '-file' => test_input_file('ecolitst.bls'));

$result = $searchio->next_result;

is($result->database_name, 'ecoli.aa', 'database_name()');
is($result->database_entries, 4289);
is($result->database_letters, 1358990);

is($result->algorithm, 'BLASTP');
like($result->algorithm_version, qr/^2\.1\.3/);
like($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
is($result->query_length, 820);
is($result->get_statistic('kappa'), '0.135');
is($result->get_statistic('kappa_gapped'), '0.0410');
is($result->get_statistic('lambda'), '0.319');
is($result->get_statistic('lambda_gapped'), '0.267');
is($result->get_statistic('entropy'), '0.383');
is($result->get_statistic('entropy_gapped'), '0.140');

is($result->get_statistic('dbletters'), 1358990);
is($result->get_statistic('dbentries'), 4289);
is($result->get_statistic('effective_hsplength'), 47);
is($result->get_statistic('effectivespace'), 894675611);
is($result->get_parameter('matrix'), 'BLOSUM62');
is($result->get_parameter('gapopen'), 11);
is($result->get_parameter('gapext'), 1);
is($result->get_statistic('S2'), '92');
is($result->get_statistic('S2_bits'), '40.0');
float_is($result->get_parameter('expect'), '1.0e-03');
is($result->get_statistic('num_extensions'), '82424');


@valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113', '0', 1567],
	      [ 'gb|AAC76922.1|', 810, 'AAC76922', '1e-91', 332],
	      [ 'gb|AAC76994.1|', 449, 'AAC76994', '3e-47', 184]);
$count = 0;
while (my $hit = $result->next_hit) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    float_is($hit->significance, shift @$d);
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while (my $hsp = $hit->next_hsp) {
            is($hsp->query->start, 1);
            is($hsp->query->end, 820);
            is($hsp->hit->start, 1);
            is($hsp->hit->end, 820);
            is($hsp->length('hsp'), 820);
            is($hsp->start('hit'), $hsp->hit->start);
            is($hsp->end('query'), $hsp->query->end);
            is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
            float_is($hsp->evalue, 0.0);
            is($hsp->score, 4058);
            is($hsp->bits,1567);	    	    
            is(sprintf("%.2f",$hsp->percent_identity), 98.29);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.9829);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.9829);
            is($hsp->gaps, 0);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

$searchio = Bio::SearchIO->new('-format' => 'blast_pull',
							  '-file'   => test_input_file('a_thaliana.blastn'));

$result = $searchio->next_result;
is($result->database_name, 'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS, GSS,or phase 0, 1 or 2 HTGS sequences)');
is($result->database_letters, 4677375331);
is($result->database_entries, 1083200);
is($result->algorithm, 'BLASTN');
like($result->algorithm_version, qr/^2\.2\.1/);
is($result->query_name, '');
is($result->query_length, 60);
is($result->get_parameter('gapopen'), 5);
is($result->get_parameter('gapext'), 2);
is($result->get_parameter('ktup'), undef);

is($result->get_statistic('lambda'), 1.37);
is($result->get_statistic('kappa'), 0.711);
is($result->get_statistic('entropy'),1.31 );
is($result->get_statistic('T'), 0);
is($result->get_statistic('A'), 30);
is($result->get_statistic('X1'), '6');
is($result->get_statistic('X1_bits'), 11.9);
is($result->get_statistic('X2'), 15);
is($result->get_statistic('X2_bits'), 29.7);
is($result->get_statistic('S1'), 12);
is($result->get_statistic('S1_bits'), 24.3);
is($result->get_statistic('S2'), 17);
is($result->get_statistic('S2_bits'), 34.2);

is($result->get_statistic('dbentries'), 1083200);

@valid = ( [ 'gb|AY052359.1|', 2826, 'AY052359', '3e-18', 96, 1, 60, 
	     '1.0000'],
	   [ 'gb|AC002329.2|AC002329', 76170, 'AC002329', '3e-18', 96, 1, 60, 
	     '1.0000' ],
	   [ 'gb|AF132318.1|AF132318', 5383, 'AF132318', '0.04', 42, 35, 55, 
	     '0.3500']);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    float_is($hit->significance, shift @$d);
    is($hit->raw_score, shift @$d );
    is($hit->start, shift @$d);
    is($hit->end,shift @$d);    
    is(sprintf("%.4f",$hit->frac_aligned_query), shift @$d);
    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1);
            is($hsp->query->end, 60);
            is($hsp->query->strand, 1);
            is($hsp->hit->start, 154);
            is($hsp->hit->end, 212);
            is($hsp->hit->strand, 1);
            is($hsp->length('hsp'), 60);	    
            float_is($hsp->evalue, 3e-18);
            is($hsp->score, 48);
            is($hsp->bits,95.6);
            is(sprintf("%.2f",$hsp->percent_identity), 96.67);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.9667);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.9831);
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->query->seq_id, '');
            is($hsp->hit->seq_id, 'gb|AY052359.1|');
            is($hsp->gaps('query'), 0);
            is($hsp->gaps('hit'), 1);
            is($hsp->gaps, 1);
            is($hsp->query_string, 'aggaatgctgtttaattggaatcgtacaatggagaatttgacggaaatagaatcaacgat');
            is($hsp->hit_string, 'aggaatgctgtttaattggaatca-acaatggagaatttgacggaaatagaatcaacgat');
            is($hsp->homology_string, '|||||||||||||||||||||||  |||||||||||||||||||||||||||||||||||');
            is(sprintf("%.2f",$hsp->get_aln->overall_percentage_identity), 96.67);
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity), 98.31);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
} 
is(@valid, 0);

$searchio = Bio::SearchIO->new('-format' => 'blast_pull',
							  '-file'   => test_input_file('frac_problems.blast'));
my @expected = ("1.000", "0.943");
while (my $result = $searchio->next_result) {
    my $hit = $result->next_hit;
	if (@expected == 2) {
		is($hit->frac_identical, shift @expected);
	}
	else {
		TODO: {
			local $TODO = 'frac_identical failing!';
			is($hit->frac_identical, shift @expected);
		}
	}
}
is(@expected, 0);

# And even more: frac_aligned_query should never be over 1!
$searchio = Bio::SearchIO->new('-format' => 'blast_pull',
							  '-file'   => test_input_file('frac_problems2.blast'));
$result = $searchio->next_result;
$hit = $result->next_hit;
is $hit->frac_aligned_query, 0.97;

# Also, start and end should be sane
$searchio = Bio::SearchIO->new('-format' => 'blast_pull',
							  '-file'   => test_input_file('frac_problems3.blast'));
$result = $searchio->next_result;
$hit = $result->next_hit;
is $hit->start('sbjct'), 207;
is $hit->end('sbjct'), 1051;

# Do a multiblast report test
$searchio = Bio::SearchIO->new('-format' => 'blast_pull',
							   '-file'   => test_input_file('multi_blast.bls'));

@expected = qw(CATH_RAT CATL_HUMAN CATL_RAT PAPA_CARPA);
my $results_left = 4;
while( my $result = $searchio->next_result ) {
    is($result->query_name, shift @expected, "Multiblast query test");
    $results_left--;
}
is($results_left, 0);

# Web blast result parsing
$searchio = Bio::SearchIO->new(-format => 'blast_pull',
							   -file   => test_input_file('catalase-webblast.BLASTP'));
ok($result = $searchio->next_result);
ok($hit = $result->next_hit);
is($hit->name, 'gi|40747822|gb|EAA66978.1|', 'full hit name');
is($hit->accession, 'EAA66978', 'hit accession');
ok(my $hsp = $hit->next_hsp);
is($hsp->query->start, 1, 'query start');
is($hsp->query->end, 528, 'query start');

# tests for new BLAST 2.2.13 output
$searchio = Bio::SearchIO->new(-format => 'blast_pull',
							  -file   => test_input_file('new_blastn.txt'));

$result = $searchio->next_result;
is($result->database_name, 'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS,GSS,environmental samples or phase 0, 1 or 2 HTGS sequences)');
is($result->database_entries, 3742891);
is($result->database_letters, 16670205594);
is($result->algorithm, 'BLASTN');
is($result->algorithm_version, '2.2.13 [Nov-27-2005]');
is($result->query_name, 'pyrR,');
is($result->query_length, 558);
is($result->get_statistic('kappa'), '0.711');
is($result->get_statistic('kappa_gapped'), '0.711');
is($result->get_statistic('lambda'), '1.37');
is($result->get_statistic('lambda_gapped'), '1.37');
is($result->get_statistic('entropy'), '1.31');
is($result->get_statistic('entropy_gapped'), '1.31');
is($result->get_statistic('dbletters'), '-509663586');
is($result->get_statistic('dbentries'), 3742891);
is($result->get_statistic('effective_hsplength'), undef);
is($result->get_statistic('effectivespace'), 8935230198384);
is($result->get_parameter('matrix'), 'blastn matrix:1 -3');
is($result->get_parameter('gapopen'), 5);
is($result->get_parameter('gapext'), 2);
is($result->get_statistic('S2'), '60');
is($result->get_statistic('S2_bits'), '119.4');
float_is($result->get_parameter('expect'), '1e-23');
is($result->get_statistic('num_extensions'), '117843');

@valid = ( [ 'gi|41400296|gb|AE016958.1|', 4829781, 'AE016958', '6e-059', 236],
	      [ 'gi|54013472|dbj|AP006618.1|', 6021225, 'AP006618', '4e-026', 127],
	      [ 'gi|57546753|dbj|BA000030.2|', 9025608, 'BA000030', '1e-023', 119]);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    float_is($hit->significance, shift @$d);
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 262);
            is($hsp->query->end, 552);
            is($hsp->hit->start, 1166897);
            is($hsp->hit->end, 1167187);
            is($hsp->length('hsp'), 291);
            is($hsp->start('hit'), $hsp->hit->start);
            is($hsp->end('query'), $hsp->query->end);
            is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
            float_is($hsp->evalue, 6e-59);
            is($hsp->score, 119);
            is($hsp->bits,236);	    	    
            is(sprintf("%.2f",$hsp->percent_identity), 85.22);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.8522);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.8522);
            is($hsp->gaps, 0);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

# Bug 2189
$searchio = Bio::SearchIO->new(-format => 'blast_pull',
							  -file   => test_input_file('blastp2215.blast'));

$result = $searchio->next_result;
is($result->database_entries, 4460989);
is($result->database_letters, 1533424333);
is($result->algorithm, 'BLASTP');
is($result->algorithm_version, '2.2.15 [Oct-15-2006]');
is($result->query_name, 'gi|15608519|ref|NP_215895.1|');
is($result->query_length, 193);
my @hits = $result->hits;
is(scalar(@hits), 10);
is($hits[1]->accession,'1W30');
float_is($hits[4]->significance,'2e-72');
is($hits[7]->score,'254');
$result = $searchio->next_result;
is($result->database_entries, 4460989);
is($result->database_letters, 1533424333);
is($result->algorithm, 'BLASTP');
is($result->algorithm_version, '2.2.15 [Oct-15-2006]');
is($result->query_name, 'gi|15595598|ref|NP_249092.1|');
is($result->query_length, 423);
@hits = $result->hits;
is(scalar(@hits), 10);
is($hits[1]->accession,'ZP_00972546');
float_is($hits[4]->significance, '0.0');
is($hits[7]->score, 624);

# Bug 2246
$searchio = Bio::SearchIO->new(-format => 'blast_pull',
                               -verbose => -1,
							  -file   => test_input_file('bug2246.blast'));
$result = $searchio->next_result;
$hit = $result->next_hit;
is $hit->name, 'UniRef50_Q9X0H5';
is $hit->length, undef;
is $hit->accession, 'UniRef50_Q9X0H5';
is $hit->description, 'Cluster: Histidyl-tRNA synthetase; n=4; Thermoto...';
is $hit->raw_score, 23;
float_is($hit->significance, 650);

#*** need to /fully/ test a multi-result, multi-hit, multi-hsp report!
