#!/usr/bin/perl

use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;
  test_begin(-tests => 7,
             -requires_modules => [qw(Statistics::Frequency
                                   Spreadsheet::ParseExcel
                                   Spreadsheet::WriteExcel)]
  );
  use_ok('Bio::Microarray::Tools::ReseqChip');
}

my $DEBUG = test_debug();

sub process_sample($$$$$$$) {

  my ($myReseqChip, $aln, $ind_id, $options_hash, $newseq_output_filename, $recalls_output_filename, $workbook) = @_;
  print "process sample $ind_id ...\n" if $DEBUG;

  $aln->sort_alphabetically;
  $myReseqChip->write_alignment2xls($aln, $workbook, $ind_id, 'human_mtDNA_RCRS', 1);
  ok(-e $workbook, 'write_alignment2xls');
  my ($newseq,$dummy);
  $newseq=$myReseqChip->calc_sequence($aln, $dummy, $options_hash, $recalls_output_filename);
  ok(-e $recalls_output_filename, 'calc_sequence');
  $myReseqChip->write2fasta($newseq, $ind_id, $newseq_output_filename, 1);
  ok(-e $newseq_output_filename, 'write2fasta');
  
}

###input files
my $Affy_frags_design_filename = test_input_file('ReseqChip_mtDNA_design_annotation_file_FINAL.xls');
my $format = 'affy_mitochip_v2';
my $Affy_sample_fasta_file = test_input_file('ReseqChip_ExampleData.fasta');
my $Mito_reference_fasta_file = test_input_file('ReseqChip_RefSeq.fasta');

#my $tmpdir = File::Spec->catfile(qw(t tmp));
#mkdir($tmpdir,0777);

###output files
my $xls_filename=test_output_file();
my $recalls_output_filename=test_output_file();
my $newseq_output_filename=test_output_file();

my $in  = Bio::SeqIO->new(-file => $Affy_sample_fasta_file,
                         -format => 'Fasta');
my $in_refseq = Bio::SeqIO->new(-file => $Mito_reference_fasta_file,
                       -format => 'Fasta');
my $refseq = $in_refseq->next_seq();
my %ref_seq_max_ins_hash=(3106 => 1);

my $myReseqChip = Bio::Microarray::Tools::ReseqChip->new($Affy_frags_design_filename, $format, \%ref_seq_max_ins_hash, $refseq);

my %options_hash=(
                include_main_sequence => 1,
                insertions => 1,
                deletions => 1,
                depth_ins => 1,
                depth_del => 9,
                depth => 1,
#                depth_ins => 0,
#                depth_del => 0,
#                depth => 0,
                consider_context => 1,
                flank_left => 10,
                flank_right => 10,
                allowed_n_in_flank => 0,
#                allowed_n_in_flank => 20,
                flank_left_ins => 4,
                flank_right_ins => 4,
                allowed_n_in_flank_ins => 1,
#                allowed_n_in_flank_ins => 21,
                flank_size_weak => 1,
                call_threshold => 55,
                ins_threshold => 35,
                del_threshold => 75,
                swap_ins => 1,
#                start_pos=> 300,
#                stop_pos=> 320,
                );


##general options
$options_hash{include_main_sequence}=1;
$options_hash{insertions}=1;
$options_hash{deletions}=1;
$options_hash{flank_size_weak}=1;
$options_hash{consider_context}=1;
$options_hash{swap_ins}=1;

##options for diploid model
$options_hash{depth}=10;
$options_hash{depth_del}=10;

$options_hash{call_threshold}=60;
$options_hash{del_threshold}=60;

$options_hash{flank_left}=1;
$options_hash{flank_right}=1;
                
$options_hash{allowed_n_in_flank}=2;

##ins options
$options_hash{depth_ins}=1;
$options_hash{ins_threshold}=60;
$options_hash{flank_left_ins}=2;
$options_hash{flank_right_ins}=2;
$options_hash{allowed_n_in_flank_ins}=2;

##options for haploid model #1
$options_hash{depth}=1;
$options_hash{depth_del}=1;

$options_hash{call_threshold}=90;
$options_hash{del_threshold}=90;

$options_hash{flank_left}=5;
$options_hash{flank_right}=5;
                
$options_hash{allowed_n_in_flank}=0;

##hap model 2
$options_hash{depth}=10;
$options_hash{depth_del}=10;

$options_hash{call_threshold}=90;
$options_hash{del_threshold}=90;

$options_hash{flank_left}=1;
$options_hash{flank_right}=1;
                
$options_hash{allowed_n_in_flank}=1;

my $ind_id="";
my $ind_id_old="";
my $workbook = Spreadsheet::WriteExcel->new($xls_filename);
my $j=1;
my $aln = Bio::SimpleAlign->new();
while ( (my $seq = $in->next_seq())) {
  my $locseq;
  my $test_complete_seq=($seq->id =~ /human_mtDNA_RCRS/);
  if ($test_complete_seq==1) {
    $ind_id=$seq->id;
  }
  else {$test_complete_seq=0;}
  if (!$test_complete_seq) {
    $locseq=$myReseqChip->insert_gaps2frag($seq);
    $aln->add_seq($locseq);
  }
  else {
    if ($aln->length>0) {
      process_sample($myReseqChip, $aln, $ind_id_old, \%options_hash,
                     $newseq_output_filename, $recalls_output_filename, $workbook);
    }
    $ind_id_old=$ind_id;
    
    $aln = new Bio::SimpleAlign();
    $locseq=$myReseqChip->insert_gaps2reference_sequence($seq);
    $aln->add_seq($locseq);
    $j++;
  }
}

process_sample($myReseqChip, $aln, $ind_id_old, \%options_hash,
               $newseq_output_filename, $recalls_output_filename, $workbook);

$workbook->close();
