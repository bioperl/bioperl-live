# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 45);
	
	use_ok('Bio::Align::DNAStatistics');
	use_ok('Bio::Align::ProteinStatistics');
	use_ok('Bio::AlignIO');
}

my $debug = test_debug();

my $in = Bio::AlignIO->new(-format => 'emboss',
			  -file   => test_input_file('insulin.water'));
my $aln = $in->next_aln();
isa_ok($aln, 'Bio::Align::AlignI');
my $stats = Bio::Align::DNAStatistics->new(-verbose => $debug);
is( $stats->transversions($aln),4);
is( $stats->transitions($aln),9);
is( $stats->pairwise_stats->number_of_gaps($aln),21);
is( $stats->pairwise_stats->number_of_comparable_bases($aln),173);
is( $stats->pairwise_stats->number_of_differences($aln),13);
is( $stats->pairwise_stats->score_nuc($aln), 224);
is( $stats->pairwise_stats->score_nuc( -aln => $aln, -match => 1,
  -mismatch => -1, -gap_open => -1, -gap_ext => -1), 126);

my $d = $stats->distance(-align => $aln,
			 -method=> 'f81');
is(  $d->get_entry('hs_insulin','seq2'), '0.07918');

$d = $stats->distance(-align=> $aln,
		      -method => 'JC');
is( $d->get_entry('hs_insulin','seq2'), '0.07918');

$d = $stats->distance(-align=> $aln,
		      -method => 'Kimura');
is( $d->get_entry('hs_insulin','seq2'), '0.07984');

$d = $stats->distance(-align=> $aln,
		      -method => 'TajimaNei');
is( $d->get_entry('seq2','hs_insulin'), '0.08106');

$d = $stats->distance(-align=> $aln,
		      -method => 'Tamura');
is( $d->get_entry('seq2','hs_insulin'), '0.08037');

#$d =  $stats->distance(-align => $aln,
#		       -method => 'JinNei');
#is( $d->get_entry('seq2','hs_insulin'), 0.0850);

$in = Bio::AlignIO->new(-format => 'clustalw',
		       -file   => test_input_file('hs_owlmonkey.aln'));

$aln = $in->next_aln();
isa_ok($aln,'Bio::Align::AlignI');

is( $stats->transversions($aln),10);
is( $stats->transitions($aln),17);
is( $stats->pairwise_stats->number_of_gaps($aln),19);
is( $stats->pairwise_stats->number_of_comparable_bases($aln),170);
is( $stats->pairwise_stats->number_of_differences($aln),27);
is( $stats->pairwise_stats->score_nuc($aln), 134);
is( $stats->pairwise_stats->score_nuc( -aln => $aln, -match => 1,
  -mismatch => -1, -gap_open => -1, -gap_ext => -1), 97);

# now test the distance calculations
$d = $stats->distance(-align => $aln, -method => 'jc');
is( $d->get_entry('human','owlmonkey'), 0.17847);

$d = $stats->distance(-align => $aln,
			  -method=> 'f81');
is(  $d->get_entry('human','owlmonkey'), '0.17847');

$d = $stats->distance(-align => $aln, -method => 'uncorrected');
is( $d->get_entry('human','owlmonkey'), 0.15882);

$d =  $stats->distance(-align => $aln, -method => 'Kimura');
is( $d->get_entry('human','owlmonkey'), 0.18105);

$d =  $stats->distance(-align => $aln, -method => 'TajimaNei');
is( $d->get_entry('human','owlmonkey'), 0.18489);

$d =  $stats->distance(-align => $aln,
			   -method => 'Tamura');

is( $d->get_entry('human','owlmonkey'), 0.18333);
#$d =  $stats->distance(-align => $aln,
#		       -method => 'JinNei');
#is( $d->get_entry('human','owlmonkey'), 0.2079);

### now test Nei_gojobori methods, hiding the expected warnings so we can
# avoid printing them ###
$stats->verbose($debug ? $debug : -1);
my ($alnobj, $result);
$in = Bio::AlignIO->new(-format => 'fasta',
			-file   => test_input_file('nei_gojobori_test.aln'));
$alnobj = $in->next_aln();
isa_ok($alnobj,'Bio::Align::AlignI');
$result = $stats->calc_KaKs_pair($alnobj, 'seq1', 'seq2');
is (sprintf ("%.1f", $result->[0]{'S'}), 40.5);
is (sprintf ("%.1f", $result->[0]{'z_score'}), '4.5');
$result = $stats->calc_all_KaKs_pairs($alnobj);
is (int( $result->[1]{'S'}), 41);
is (int( $result->[1]{'z_score'}), 4);
$result = $stats->calc_average_KaKs($alnobj, 100);
is (sprintf ("%.4f", $result->{'D_n'}), 0.1628);
$stats->verbose($debug);

# now test Protein Distances
my $pstats = Bio::Align::ProteinStatistics->new();
$in = Bio::AlignIO->new(-format => 'clustalw',
			-file   => test_input_file('testaln.clustalw'));
$alnobj = $in->next_aln();
isa_ok($alnobj,'Bio::Align::AlignI');
$result = $pstats->distance(-method => 'Kimura',
			    -align  => $alnobj);
isa_ok($result, 'Bio::Matrix::PhylipDist');

is ($result->get_entry('P84139','P814153'),   '0.01443');
is ($result->get_entry('P841414','P851414'),  '0.01686');
is ($result->get_entry('P84139','P851414'),   '3.58352');

my $seq = Bio::Seq->new(-id=>'NOT3MUL', -seq=>'gatac');
isa_ok($seq, 'Bio::PrimarySeqI');
eval { 
  Bio::Align::DNAStatistics->count_syn_sites($seq); 
};
like($@, qr/not integral number of codons/);

# bug 2901
$in = Bio::AlignIO->new(-file => test_input_file('bug2901.fa'),
                        -format => 'fasta');

$stats = Bio::Align::DNAStatistics->new(-verbose => 2);
$aln = $in->next_aln();
my $matrix;
throws_ok {
$matrix = $stats->distance(-align=>$aln,-method=>'Uncorrected');
} qr/No distance calculated between seq3 and seq4/, "Warn if seqs don't overlap";
$stats->verbose(-1);
$matrix = $stats->distance(-align=>$aln,-method=>'Uncorrected');
like($matrix->print_matrix, qr/-1/);

