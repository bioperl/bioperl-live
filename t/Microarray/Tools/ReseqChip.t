#!/usr/bin/perl

use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;
  test_begin(-tests => 10,
             -requires_modules => [qw(Statistics::Frequency
                                   Spreadsheet::ParseExcel
                                   Spreadsheet::WriteExcel)]
  );
  use_ok('Bio::Microarray::Tools::ReseqChip');
}

my $DEBUG = test_debug();

sub read_params($$$$$) {

  my ($params_filename, $options_hash, $ncall_threshold, $ploidy, $mut_type) = @_;
  
  open(PARAMSFILE, $params_filename);
  if (!(-e $params_filename)) {
    print "File ($params_filename) does not exists !!!\n";
  }
  my $old_ncall=0;
  my $cur_ncall=0;
  my @params;
  my @p_array;
  my $test=0;
  while (<PARAMSFILE>) {
    if (/$ploidy/) {
      if (/$mut_type/) {
        $test=1;
        @params = split(' ', $_);
        @p_array=split('_', $params[0]);
        $cur_ncall=$params[1];
        if (($ncall_threshold<$cur_ncall and $ncall_threshold>=$old_ncall)) {
          last;
        }
        $old_ncall=$cur_ncall;
      }
    }
  }
  close(PARAMSFILE);
  #print("the ncalllevel: $cur_ncall\n");
  ok($test, 'read_params');

  if ($mut_type eq 'Subdel') {
    $options_hash->{depth}=$p_array[3];
    $options_hash->{depth_del}=$p_array[3];
    $options_hash->{flank_left}=$p_array[4];
    $options_hash->{flank_right}=$p_array[5];
    $options_hash->{allowed_n_in_flank}=$p_array[6];
    $options_hash->{call_threshold}=$p_array[7];
    $options_hash->{del_threshold}=$p_array[7];
  } else {
    $options_hash->{depth_ins}=$p_array[3];
    $options_hash->{flank_left_ins}=$p_array[4];
    $options_hash->{flank_right_ins}=$p_array[5];
    $options_hash->{allowed_n_in_flank_ins}=$p_array[6];
    $options_hash->{ins_threshold}=$p_array[7];  
  }
    
}


sub process_sample($$$$$$$) {

  my ($myReseqChip, $aln, $ind_id, $options_hash, $newseq_output_filename, $recalls_output_filename, $workbook) = @_;
  print "process sample $ind_id ...\n" if $DEBUG;

  $aln->sort_alphabetically;
  $myReseqChip->write_alignment2xls($aln, $workbook, $ind_id, 'human_mtDNA_RCRS', 1);

  my $newseq=$myReseqChip->calc_sequence($aln, $options_hash, $recalls_output_filename);
  ##test if $newseq has expected length
  ok(length($newseq)==16565, 'calc_sequence: length');
  ##test if logfile has expected size
  #print((-s $recalls_output_filename)."logfile size\n");
  ok(((-s $recalls_output_filename)==9242 or (-s $recalls_output_filename)==20622), 'calc_sequence: logfile');

  $myReseqChip->write2fasta($newseq, $ind_id, $newseq_output_filename, 1);
  ##test if fastafile has expected size
  #print((-s $newseq_output_filename)."fastafile size\n");
  ok(((-s $newseq_output_filename)==16873 or (-s $newseq_output_filename)==33746), 'write2fasta');
  
}

###input files
my $Affy_frags_design_filename = test_input_file('ReseqChip_mtDNA_design_annotation_file_FINAL.xls');
my $format = 'affy_mitochip_v2';
my $Affy_sample_fasta_file = test_input_file('ReseqChip_ExampleData.fasta');
my $Affy_alternative_sample_fasta_file = test_input_file('ReseqChip_ExampleData.fasta');
my $Mito_reference_fasta_file = test_input_file('ReseqChip_RefSeq.fasta');
my $Parameter_file = test_input_file('ReseqChip_ParamsNcall.csv');

#my $tmpdir = File::Spec->catfile(qw(t tmp));
#mkdir($tmpdir,0777);

###output files
my $xls_filename=test_output_file();
my $recalls_output_filename=test_output_file();
my $newseq_output_filename=test_output_file();

my $in = Bio::SeqIO->new(-file => $Affy_sample_fasta_file,
                         -format => 'Fasta');

my $in_alt  = Bio::SeqIO->new(-file => $Affy_alternative_sample_fasta_file,
                         -format => 'Fasta');

my $in_refseq = Bio::SeqIO->new(-file => $Mito_reference_fasta_file,
                       -format => 'Fasta');
my $refseq = $in_refseq->next_seq();
my %ref_seq_max_ins_hash=(3106 => 1);

my $myReseqChip = Bio::Microarray::Tools::ReseqChip->new($Affy_frags_design_filename, $format, \%ref_seq_max_ins_hash, $refseq);


##general options
my %options_hash=(
  include_main_sequence => 1,
  insertions => 1,
  deletions => 1,
  swap_ins => 1,
  flank_size_weak => 1,
  consider_context => 1
);


##data specific options have to set by parsing parameter file

##subdel
my $subdel_ncalls=4.5;  #specify value between 0 and 16.9 for Hap
                        # and 0 and 55.6 for Dip model, respectively
read_params($Parameter_file, \%options_hash, $subdel_ncalls, 'Hap', 'Subdel');

##insertions
my $ins_ncalls=5.4; #specify value between 0 and 20.9
                    # and 0.1 and 54.9
read_params($Parameter_file, \%options_hash, $ins_ncalls, 'Hap', 'Ins');
#for my $pos (sort{$a<=>$b}keys %options_hash) {
#  print "$pos :".$options_hash{$pos}."\n";
#}
            


my $ind_id="";
my $ind_id_old="";
my $workbook = Spreadsheet::WriteExcel->new($xls_filename);
my $j=1;
my $aln = Bio::SimpleAlign->new();
while ( (my $seq = $in->next_seq())) {
  my ($locseq, $locseq_alt);
  my $seq_alt=$in_alt->next_seq();
  my $test_complete_seq=($seq->id =~ /human_mtDNA_RCRS/);
  if ($test_complete_seq==1) {
    $ind_id=$seq->id;
  }
  else {$test_complete_seq=0;}
  if (!$test_complete_seq) {
    $locseq=$myReseqChip->insert_gaps2frag($seq);
    $aln->add_seq($locseq);
    
    ##alternative basecalls
    $locseq_alt=$myReseqChip->insert_gaps2frag($seq_alt);
    $options_hash{alternative_sequence_hash}->{$locseq_alt->id}=$locseq_alt;
    
  } else {
    if ($aln->length>0) {
      process_sample($myReseqChip, $aln, $ind_id_old, \%options_hash,
                     $newseq_output_filename, $recalls_output_filename, $workbook);
    }
    $ind_id_old=$ind_id;
    
    $aln = new Bio::SimpleAlign();
    $locseq=$myReseqChip->insert_gaps2reference_sequence($seq);
    $aln->add_seq($locseq);

    ##alternative primary basecalls for insertions
    $locseq_alt=$myReseqChip->insert_gaps2reference_sequence($seq_alt);
    $options_hash{alternative_sequence_hash}->{$locseq_alt->id}=$locseq_alt;
    
    $j++;
  }
}

process_sample($myReseqChip, $aln, $ind_id_old, \%options_hash,
               $newseq_output_filename, $recalls_output_filename, $workbook);

$workbook->close();
##test if xls file has expected size
#print((-s $xls_filename)."xlsfile size\n");
cmp_ok((-s $xls_filename), '>', 1100000, 'write_alignment2xls');
