# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

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
    plan tests => 36; 
}

if( $error == 1 ) {
    exit(0);
}

my $debug = -1;

use Bio::Align::DNAStatistics;
use Bio::Align::ProteinStatistics;
use Bio::AlignIO;
use Bio::Root::IO;

my $in = new Bio::AlignIO(-format => 'emboss',
			  -file   => Bio::Root::IO->catfile('t', 'data',
							    'insulin.water'));
my $aln = $in->next_aln();
ok($aln);
my $stats = new Bio::Align::DNAStatistics(-verbose => $debug);
ok( $stats->transversions($aln),4);
ok( $stats->transitions($aln),9);
ok( $stats->pairwise_stats->number_of_gaps($aln),21);
ok( $stats->pairwise_stats->number_of_comparable_bases($aln),173);
ok( $stats->pairwise_stats->number_of_differences($aln),13);

my $d = $stats->distance(-align => $aln,
			 -method=> 'f81');
ok(  $d->get_entry('hs_insulin','seq2'), '0.07918');

$d = $stats->distance(-align=> $aln,
		      -method => 'JC');
ok( $d->get_entry('hs_insulin','seq2'), '0.07918');

$d = $stats->distance(-align=> $aln,
		      -method => 'Kimura');
ok( $d->get_entry('hs_insulin','seq2'), '0.07984');

$d = $stats->distance(-align=> $aln,
		      -method => 'TajimaNei');
ok( $d->get_entry('seq2','hs_insulin'), '0.08106');

$d = $stats->distance(-align=> $aln,
		      -method => 'Tamura');
ok( $d->get_entry('seq2','hs_insulin'), '0.08037');

#$d =  $stats->distance(-align => $aln,
#		       -method => 'JinNei');
#ok( $d->get_entry('seq2','hs_insulin'), 0.0850);

$in = new Bio::AlignIO(-format => 'clustalw',
		       -file   => Bio::Root::IO->catfile('t','data',
							 'hs_owlmonkey.aln'));

$aln = $in->next_aln();
ok($aln);

ok( $stats->transversions($aln),10);
ok( $stats->transitions($aln),17);
ok( $stats->pairwise_stats->number_of_gaps($aln),19);
ok( $stats->pairwise_stats->number_of_comparable_bases($aln),170);
ok( $stats->pairwise_stats->number_of_differences($aln),27);

# now test the distance calculations
$d = $stats->distance(-align => $aln, -method => 'jc');
ok( $d->get_entry('human','owlmonkey'), 0.17847);

$d = $stats->distance(-align => $aln,
		      -method=> 'f81');
ok(  $d->get_entry('human','owlmonkey'), '0.17847');

$d = $stats->distance(-align => $aln, -method => 'uncorrected');
ok( $d->get_entry('human','owlmonkey'), 0.15882);

$d =  $stats->distance(-align => $aln, -method => 'Kimura');
ok( $d->get_entry('human','owlmonkey'), 0.18105);

$d =  $stats->distance(-align => $aln, -method => 'TajimaNei');
ok( $d->get_entry('human','owlmonkey'), 0.18489);

$d =  $stats->distance(-align => $aln,
		       -method => 'Tamura');

ok( $d->get_entry('human','owlmonkey'), 0.18333);
#$d =  $stats->distance(-align => $aln,
#		       -method => 'JinNei');
#ok( $d->get_entry('human','owlmonkey'), 0.2079);

### now test Nei_gojobori methods ##
$in = Bio::AlignIO->new(-format => 'fasta',
			-file   => Bio::Root::IO->catfile('t','data',
							  'nei_gojobori_test.aln'));
my $alnobj = $in->next_aln();
ok($alnobj);
my $result = $stats->calc_KaKs_pair($alnobj, 'seq1', 'seq2');
ok (sprintf ("%.1f", $result->[0]{'S'}), 40.5);
ok (sprintf ("%.1f", $result->[0]{'z_score'}), '4.5');
$result = $stats->calc_all_KaKs_pairs($alnobj);
ok (int( $result->[1]{'S'}), 41);
ok (int( $result->[1]{'z_score'}), 4);
$result = $stats->calc_average_KaKs($alnobj, 100);
ok (sprintf ("%.4f", $result->{'D_n'}), 0.1628);


# now test Protein Distances
my $pstats = Bio::Align::ProteinStatistics->new();
$in = Bio::AlignIO->new(-format => 'clustalw',
			-file   => Bio::Root::IO->catfile('t','data',
							  'testaln.aln'));
$alnobj = $in->next_aln();
ok($alnobj);
$result = $pstats->distance(-method => 'Kimura',
			    -align  => $alnobj);
ok($result);

ok ($result->get_entry('P84139','P814153'),   '0.01443');
ok ($result->get_entry('P841414','P851414'),  '0.01686');
ok ($result->get_entry('P84139','P851414'),   '3.58352');

my $seq = Bio::Seq->new(-id=>'NOT3MUL', -seq=>'gatac');
ok($seq);
eval { 
  Bio::Align::DNAStatistics->count_syn_sites($seq); 
};
ok($@ =~ m/not integral number of codons/);


