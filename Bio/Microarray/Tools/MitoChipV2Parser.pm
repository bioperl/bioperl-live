#------------------------------------------------------------------------------
# PACKAGE : Bio::Microarray::Tools::MitoChipV2Parser.pm
# PURPOSE : Parser for affy design file
# AUTHOR  : Marian Thieme
# CREATED : 28.09.2007
# REVISION: 
# STATUS  : Beta
#
# For documentation, run this module through pod2html
#------------------------------------------------------------------------------


=head1 NAME

Bio::Microarray::Tools::MitoChipV2Parser- Class for parsing design file for Affy
MitoChip V2.0

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 AUTHORS
        
Marian Thieme (marian.thieme@arcor.de)

=head1 COPYRIGHT
        
Copyright (c) 2007 Institute of Functional Genomics, University Regensburg,
granted by Baygene. All Rights Reserved.
        
This module is free software; you can redistribute it and/or modify it under the
same terms as Perl itself.

=head1 DISCLAIMER
        
This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::Microarray::Tools::MitoChipV2Parser;

use Spreadsheet::ParseExcel;                                   

use strict;
use warnings;

use base qw(Bio::Root::Root);



=head2 new()
  Title		: new/Creator method
  Usage		: my $parser=Bio::Microarray::Tools::MitoChipV2Parser->new($file_name);
  Function	: Creates 2 hashes by parsing the specified design file. 
  Returns	: MitoChipV2Parser object, consists mainly of 2 hashes:
                  1. $self->{_frags_hash} hash of hashs containing the fragment id as key and a hash as value. 
                  The nested hash has a position (covering a variying site hold by the fragment) as key and an 
                  array of length 4 (a,b,c,d) containing:
                    a. the mutation type (ins, del, sub)
                    b. the mutated base (A, C, G, T (in the case of subsititution or insertion of length 1), 
                        D (in the case of a deletion) 
                        and for instance CC, CCC, ... in the case of an insertion of length 2 or 3)
                    c. the start position of the fragment
                    d. the end position of the fragment (that information is identical for each array belonging 
                    to the same fragment)
                  2. $self->{_oligos2calc_hash} hash containing the start and stop position of location of fragments 
                  (key: start pos, value: end pos)
  Args		: $filename - name of affy mitochip v2 design file (for instance: mtDNA_design_annotation_file_FINAL.xls)

=cut


sub new {

  my ($class, $file_name) = @_;
  my $self = $class->SUPER::new();
  $self->{_frags_hash}=undef;
  $self->{_oligos2calc_hash}=undef;
  $self->throw("Must provide filename as first argument !") unless $file_name;

  
  my %max_ins_hash=();
  $self->{_frags_hash}=$self->_parse_Affy_mtDNA_design_annotation_file($file_name);
  #print "Created array of redundant fragments.\n";
  #print $self->{_frags_hash};
  bless $self, $class;
  return $self;
}

=head2 _parse_Affy_mtDNA_design_annotation_file()

 Title     : _parse_Affy_mtDNA_design_annotation_file()
 Usage     : private Function, dont call directly
 Function  : Creates Hash of insertions of maximal length, by parsing the affy chip 
             design file and calcs the location of clusters of the Redundant Fragments
             
 Returns   : Returns hash of maximal insertions and calls private function _calc_oligo_region_hash()
 Args      : $filename - excel filename of the affy design file

=cut

sub _parse_Affy_mtDNA_design_annotation_file() {

  my ($self, $filename) = @_;
  my $oExcel = new Spreadsheet::ParseExcel;
  my $oBook = $oExcel->Parse($filename);
  $self->throw("File $filename not found") unless $oBook;
  #my ($iR, $iC, $oWkS, $oWkC);
  
  my %frags=();
  #309.1C_309.2C_309.3C_315.1C_337D is parsed to hash:
  #mt228ref285-339 => 	315 =>  ins ... C ... 285 .. 339 .
  #			309 =>  del ... D ... 285 .. 339
  #mt211ref266-339 => 	315 =>  ins ... C ... 266 .. 339 .
  #			309 =>  ins ... CC ... 266 .. 339 .
  #			291 =>  del ... D ... 266 .. 339 .
  #			290 =>  del ... D ... 266 .. 339 .
  my %oligo_region_hash=(0 => 0);
  my %hhash;
  my $oWkS = $oBook->{Worksheet}[0];
  
  my $count_row=0;
  my $count_var=0;
  my $count_all_var=0;
  for(my $iR = $oWkS->{MinRow}+1 ; defined $oWkS->{MaxRow} && ($iR <= $oWkS->{MaxRow}) ; $iR++)	{
    my $fragment_id = $oWkS->{Cells}[$iR][$oWkS->{MinCol}];
    my $mutation_desc = $oWkS->{Cells}[$iR][$oWkS->{MinCol}+1];
    my $frag_pos = $oWkS->{Cells}[$iR][$oWkS->{MinCol}+2];
    my $mito_pos = $oWkS->{Cells}[$iR][$oWkS->{MinCol}+3];
    my $insertion= $oWkS->{Cells}[$iR][$oWkS->{MinCol}+4];
    my @array;
    #print $fragment_id->Value."\n";
    @array= ($fragment_id->Value =~ m/(\d+)/g);
    my $start=$array[1];
    my $stop=$array[2];
    
    %oligo_region_hash=$self->_calc_oligo_region_hash($start, $stop, \%oligo_region_hash);
    
    
    my @mut_array=split("_",$mutation_desc->Value);
    #print "@mut_array\n";
    push(@mut_array,'dummy');
    my @pos_array;
    my @mut_pos=($mut_array[0] =~ m/(\d+)/g);
    my $mut_pos=$mut_pos[0];
    my @res;
    my $i=0;
    my $length=@mut_array;
    #foreach my $item (@mut_array) {
    #  print "$item\n";
    #}
    foreach my $item (@mut_array) {
      my $test_write=1;
      $i++;
      if ($item eq 'dummy') {
        @array=(0);
      }
      else {
        @array= ($item =~ m/(\d+)/g);
      }
      
      #print @array." @array\n";
      if (exists($frags{$fragment_id->Value}{$array[0]})) {
        $res[1]=$res[1].substr($item,-1);
      }
      else {
        if ($array[0] != $mut_pos or $i==$length) {
          $frags{$fragment_id->Value}{$mut_pos}=[@res];
          $test_write=0;
          @res=();
          $count_var++;
        }
        if ($i!=$length) {
          if ($item =~ /\./) { @res=('ins');}
          elsif (substr($item,-1) eq "D") {@res=('del')}
          else { @res=('sub');}
          push (@res, substr($item,-1), $start, $stop);
          $frags{$fragment_id->Value}{$array[0]}=();
        }
      }
      
      $mut_pos=$array[0];
    }
    $count_row++;
    $count_all_var+=$i;
    #print "Total no of add. variants for current fragment ($start $stop): $i ($count_all_var)\n";
  }
  #print "$count_row fragments having $count_var variants\n";
  $self->{_oligos2calc_hash}=$self->_check_overlapping_regions(\%oligo_region_hash);
  #$self->count_different_variants(\%frags);
  return \%frags;
}


=head2 count_different_variants()

 Title     : count_different_variants()
 Usage     : MyMitochipParser->count_different_variants($myfrags) (For Developder only)
 Function  : Counts number of Positions covered by alternative probes, number of total and single (unique) variants
             separately for single pointmutations, deletions and insertions.
 Args      : Reference of Hash of alternative Probes (which is provided by _parse_Affy_mtDNA_design_annotation_file())

=cut

sub count_different_variants() {
  
  my ($self, $frags) = @_;
  my %diff_var_hash=();
  foreach my $key (sort keys %{$frags}) {
    #print $key.": ";
    #print $frags->{$key};
    foreach my $subkey (sort keys %{$frags->{$key}}) {
      #print "\n$subkey:";
      my @array=@{$frags->{$key}{$subkey}};
      #print $array[0]." ".$array[1];
      #print "\n";
      # $type=sub, del, ins
      # $var= A, CCC, --, ...
      my $pos=$subkey;
      my $type=$array[0];
      my $var=$array[1];
      if (exists($diff_var_hash{$pos})) {
        if (exists($diff_var_hash{$pos}{$type})) {
          if (exists($diff_var_hash{$pos}{$type}{$var})) {
            $diff_var_hash{$pos}{$type}{$var}++;
          }
          else {
            $diff_var_hash{$pos}{$type}{$var}=1;
          }
        }
        else {
          $diff_var_hash{$pos}{$type}{$var}=1;
        }
        
      }
      #insert position
      else {
        #print "start: $pos $type $var =1\n";
        $diff_var_hash{$pos}{$type}{$var}=1;
      }
      
    }
    
    #print "\n\n";
    
  }
  my $count_pos=0;
  my $count_type=0;
  my $count_var=0;
  my $count_var_ins=0;
  my $count_var_del=0;
  my $count_var_sub=0;
  my $count_tot=0;
  my $count_var_ins_tot=0;
  my $count_var_del_tot=0;
  my $count_var_sub_tot=0;
  my $no=0;
  foreach my $pos (sort{$a<=>$b} keys %diff_var_hash) {
    $no++;
    $count_pos++;
    print "\n$no: $pos:";
    foreach my $type (keys %{$diff_var_hash{$pos}}) {

      $count_type++;
      print " $type:";
      foreach my $var (keys %{$diff_var_hash{$pos}{$type}}) {

        $count_var++;
        $count_tot+=$diff_var_hash{$pos}{$type}{$var};
        if ($type eq "ins") {
          $count_var_ins++;
          $count_var_ins_tot+=$diff_var_hash{$pos}{$type}{$var};
        }
        elsif ($type eq "del") {
          $count_var_del++;
          $count_var_del_tot+=$diff_var_hash{$pos}{$type}{$var};
        }
        elsif ($type eq "sub") {
          $count_var_sub++;
          $count_var_sub_tot+=$diff_var_hash{$pos}{$type}{$var};
        
        }
        else {print "doesnt exist!\n";}
        print " $var: $count_var ($count_tot)";
      }
    }
  }
  print "\n\n$count_pos/$count_type/$count_var/$count_tot\n\n";
  print "ins: $count_var_ins $count_var_ins_tot\n";
  print "del: $count_var_del $count_var_del_tot\n";
  print "sub: $count_var_sub $count_var_sub_tot\n";
}

=head2 _calc_oligo_region_hash()

 Title     : _calc_oligo_region_hash()
 Usage     : private Function, dont call directly
 Function  : Adds redundant fragment to hash of location of clusters of the redundant fragments
 Returns   : oligo_region_hash, calculated so far
 Args      : $start - start position of current fragment
             $stop - stop position of current fragment
             $oligo_region_hash (reference) - hash of cluster of oligo regions
=cut

sub _calc_oligo_region_hash() {
  
  my ($self, $start, $stop, $oligo_region_hash) = @_;
    
    my %hhash=%{$oligo_region_hash};
    if (exists($oligo_region_hash->{0})) {
      delete($oligo_region_hash->{0});
      $oligo_region_hash->{$start}=$stop;
    }
    else {
      my $test_new=1;
      foreach my $key (sort{$a<=>$b} keys %{$oligo_region_hash}) {
        if ($oligo_region_hash->{$key}>=$stop and $key<=$start) {
          $test_new=0;
        }
        elsif ($key <= $start and $oligo_region_hash->{$key}>= $start and $stop-1 >= $oligo_region_hash->{$key}) {
          $hhash{$key}=$stop;
          $test_new=0;
        }
        elsif ($oligo_region_hash->{$key}>=$stop and $key<=$stop and $start+1<=$key) {
          delete($hhash{$key});
          $hhash{$start}=$oligo_region_hash->{$key};
          $test_new=0;          
        }
        ##both ends are outside the current
        elsif ($oligo_region_hash->{$key}<=$stop and $key>=$start) {
          $test_new=0;
        }
      }
      if ($test_new) {
        $oligo_region_hash->{$start}=$stop;        
      }
      else {
        $oligo_region_hash=\%hhash;
      }
    }
    return %{$oligo_region_hash};
}

=head2 _check_overlapping_regions()

 Title     : _check_overlapping_regions()
 Usage     : private Function, dont call directly
 Function  : checks for overerlapping clusters Redundant Fragments and combines
             they to one if those are found
 Returns   : (cleared) oligo_region_hash
 Args      : $oligo_region_hash (reference) - hash of cluster of oligo regions


=cut

sub _check_overlapping_regions() {
  my ($self, $oligo_region_hash) = @_;
  
  my %hhash=%{$oligo_region_hash};
  my $oldkey=0;
  my $oldvalue=0;
  foreach my $key (sort{$a<=>$b} keys %{$oligo_region_hash}) {
    if ($oldkey+1<=$key and $oldvalue>=$key and $oligo_region_hash->{$key}>=$oldvalue) {
      delete($hhash{$key});
      $hhash{$oldkey}=$oligo_region_hash->{$key};
    }
    $oldkey=$key;
    $oldvalue=$oligo_region_hash->{$key};
  }
  return \%hhash;
}


1;