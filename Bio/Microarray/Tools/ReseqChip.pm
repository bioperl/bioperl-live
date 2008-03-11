#------------------------------------------------------------------------------
# PACKAGE : Bio::Microarray::Tools::ReseqChip
# PURPOSE : Analyse redundant fragments of Affymetrix Resequencing Chip
# AUTHOR  : Marian Thieme
# CREATED : 21.09.2007
# REVISION: 
# STATUS  : Beta
#
# For documentation, run this module through pod2html
#------------------------------------------------------------------------------


=head1 NAME

Bio::Microarray::Tools::ReseqChip - Class for extraction and incorporation of information 
                                   about redundant fragments of Affy Mitochip v2.0

=head1 SYNOPSIS

 use RedundantFragments;


 my %ref_seq_max_ins_hash=(3106 => 1);
 my $reseqfragSample=Bio::Tools::ReseqChipRedundantFragments->new(
    $Affy_frags_design_filename,
    $format,
    \%ref_seq_max_ins_hash,
    $refseq_chip);

 my $aln = new Bio::SimpleAlign();
 my $in  = Bio::SeqIO->new(-file => $Affy_reseq_sample_fasta_file,
                          -format => 'Fasta');
                         
 while ( (my $seq = $in->next_seq())) {

  my $locseq;
  my $test_complete_seq=($seq->id =~ /human_mtDNA_RCRS/);
  if ($test_complete_seq==1) {
    $locseq=$reseqfragSample->insert_gaps2reference_sequence($seq);
  }
  else {
    $locseq=$reseqfragSample->insert_gaps2frag($seq);
  }
  $aln->add_seq($locseq);
 }
 my $new_sequence=$reseqfragSample->calc_sequence($aln, $options_hash [,"output_file"]);


=head1 DESCRIPTION

Process Affy MitoChip v2 Data to create an alignment of the "redundant" fragments to the reference sequence, 
taking account for insertions/deletion which are defined by Affy mtDNA_Design_Annotion.xls file. Based on the
that alignment substitutions, deletions and insertion can be detected and initally not called bases can called 
as well falsly called bases can recalled. Depending on the depth at a certain position in the alignment and 
sequence reliability (in terms of certain number of allowed Ns in a k-base-window within each redundant fragment,
contributing to a certain alignment position) initially made N Calls can called and possibly falsly called bases
can recalled. Moreover insertion and deletion as well as snps lying in highly variable regions can be detected.

Assumption:
Insertions refer to refseq, when regarding the max insertions of refseq (refseq_max_ins_hash).
Optionshash is given when calculating the sequence respect to redundant fragments

This module depends on the following modules:
use Bio::Microarray::Tools::MitoChipV2Parser
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::LiveSeq::Mutation;
use Statistics::Frequency;
use Spreadsheet::ParseExcel;
use Spreadsheet::WriteExcel;


=head1 AUTHORS
        
Marian Thieme (marian.thieme@arcor.de)

=head1 COPYRIGHT
        
Copyright (c) 2007 Institute of Functional Genomics, University Regensburg, granted by Baygene. All Rights Reserved.
        
This module is free software; you can redistribute it and/or modify it under the same terms as Perl itself.
        




=head1 DISCLAIMER
        
This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::Microarray::Tools::ReseqChip;
                                   
use vars qw(@ISA);
use strict;
use warnings;

use Bio::Root::Root;
@ISA = qw(Bio::Root::Root);

use Bio::Microarray::Tools::MitoChipV2Parser;

use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::LiveSeq::Mutation;
use Statistics::Frequency;
use Spreadsheet::ParseExcel;
use Spreadsheet::WriteExcel;


=head2 new()

 Title     : new

 Usage     : my $reseqfragSample=Bio::Resequencing::RedundantFragments->new
             ($Affy_frags_design_filename, $format, $refseq, \%oligos2calc_hash)


 Function  : Creates Hash of insertions of maximal length, by parsing the affy chip 
             design file and calcs the location of clusters of the Redundant Fragments
             by calling _parse_Affy_mtDNA_design_annotation_file() and sets further 
             member variables.
             
             
 Returns   : Returns a new RedundantFragments object
 
 Args      : $Affy_frags_design_filename (Affymetrix xls design file, 
             for instance: mtDNA_design_annotation_FINAL.xls for mitochondrial Genome)
             
             $format (only xls is available, because its that format delvered by Affymetrix)
             
             [$reseq_max_ins_hash, $refseq] (insertions as hash (pos1 => insertions length1, pos2 => insertions length2, ...) 
             for reference sequence $refseq (Locatable Sequence Object))
            
 Membervars: frags_hash		- contains all variations described by the affy_design_annotation file 
                                  (fragment_id => (pos1 => [muttype, mut, start, stop]), ... )
             max_ins_hash	- contains all (maximal) insertion, which are needed to build "alignable" 
                                  fragments by inserting appropriate gaps
             oligos2calc_hash	- describes the coverage of the redundand fragments (start => stop)
             reseq		- Locatable Sequence holds reference sequence
             refseq_max_ins_hash- insertions need to be applied to the reference sequence
             gapped_refseq	- reference sequence when insertions are done


=cut


sub new {

  my ($class, $file_name, $format, $refseq_max_ins_hash, $refseq) = @_;
  my $self = $class->SUPER::new();
  $self->{_frags_hash}=undef;
  $self->{_max_ins_hash}=undef;
  $self->{_refseq}=undef;
  $self->{_gapped_refseq}=undef;
  $self->{_oligos2calc_hash}=undef;
  $self->{_refseq_max_ins_hash}=undef;
  $self->{_frags_max_ins_hash}=undef;
  
  $self->throw("Must provide filename as first argument !") unless $file_name;

  $self->throw("Must specify format (only xls at present) as second argument !") unless $format;
  
  $self->warn("No reference sequence given !") unless $refseq;
  
  my %max_ins_hash=();
  
  
  if ($format eq "affy_mitochip_v2") {
    my $parser=Bio::Microarray::Tools::MitoChipV2Parser->new($file_name);
    #$self->{_frags_hash}=$self->_parse_Affy_mtDNA_design_annotation_file($file_name);
    $self->{_frags_hash}=$parser->{_frags_hash};
    $self->{_oligos2calc_hash}=$parser->{_oligos2calc_hash};
    print "Created array of redundant fragments. \n";
  }
  else {
    $self->warn("$format isnt supported/implemented");
  }
  
  
  if ($refseq_max_ins_hash) {
    $self->{_refseq_max_ins_hash}=$refseq_max_ins_hash;
  }
  if ($refseq) {
    $self->{_refseq}=$refseq;
    $self->{_gapped_refseq}=$self->insert_gaps2reference_sequence($refseq);
  }

  $self->{_max_ins_hash}=();
  my $correction=0;
  while ( my ($key, $value) = each %{$self->{_frags_hash}} ) {
    #print "$key: ";
    while ( my ($subkey, $subvalue) = each %$value ) {
      my @array= @{$self->{_frags_hash}->{$key}{$subkey}};
      #print "$subkey @array \n";
      if (!(exists($self->{_max_ins_hash}{$subkey+$correction}))) {
        if (@$subvalue[0] eq "ins") {
          $self->{_max_ins_hash}{$subkey+$correction}=length(@$subvalue[1]);
        }
      }	elsif ( $self->{_max_ins_hash}{$subkey+$correction} < length(@$subvalue[1]) ) {
          $self->{_max_ins_hash}{$subkey+$correction}=length(@$subvalue[1]);
      }
    }
  }
  #$self->{_max_ins_length}=$self->_get_max_ins_length();
  #print "\nMaximal insertions (Position => Length):\n";
  #foreach my $key (sort{$a<=>$b} keys %{$self->{_max_ins_hash}}) {
  #  print "$key => ".$self->{_max_ins_hash}{$key}."\n";
  #}
  my %temp_hash=%{$self->{_max_ins_hash}};
  $self->{_frags_max_ins_hash}=\%temp_hash;
  #print "\nRegions covered by redundant fragments (start position => end position):\n";
  #foreach my $key (sort{$a<=>$b} keys %{$self->{_oligos2calc_hash}}) {
  #  print "$key => ".$self->{_oligos2calc_hash}{$key}."\n";
  #}
  bless $self, $class;
  return $self;
}


=head2 _get_cont_no()

 Title     : _get_cont_no()
 Usage     : private Function, dont call directly
 Function  : returns first contiguous number from given string and cut that part away from the string
 Returns   : first contiguous number wihtin the string
 Args      : \$string : reference to a string


=cut

sub _get_cont_no() {

  my ($self, $myvalue1) = @_;
  my $myvalue=$$myvalue1;
  my $start="";
  my $first=1;
  my $thisp=1;
  while ($myvalue =~ /[0-9]/g and $thisp==1) {
    $thisp=pos($myvalue);
    if ($first == 1 or $thisp == 1) {
      $thisp=1;
      $start.= substr($myvalue,pos($myvalue)-1,1);
      $myvalue=substr($myvalue,pos($myvalue));
    }
    $first=0;
  }
  #$start.=substr($myvalue,0,1);
  if (length($myvalue)>0) {
    $myvalue=substr($myvalue,1);
    #break();
    $$myvalue1=$myvalue;
  }
  return $start;
}



=head2 _get_start_pos()

 Title     : _get_start_pos()
 Usage     : private Function, dont call directly
 Function  : calcs cumulative offset for a given position
 Returns   : position
 Args      : $pos : current position in the sequence/fragment


=cut

sub _get_start_pos() {

  my ($self, $pos, $hash, $test) = @_;
  #my $pos = (@_[1]);
  my $hashname='_max_ins_hash';
  if ($hash) {
    $hashname=$hash;
  }
  if ($test) {
    print "inside get_start_pos\n";
  }
  my $offset=0; 
  foreach my $key (sort{$a<=>$b} keys %{$self->{$hashname}}) {
    if ($test) {
      print "$key: ".$self->{$hashname}{$key}."\n";
    }
    if ($pos>$key-12) {
    #if ($pos>=$key) {
      $offset+=$self->{$hashname}{$key};
      if ($test) {
        print "($offset)\n";
      }
    }
  }
  return $offset;
}

=head2 insert_gaps2frag()

 Title     : insert_gaps2frag()
 Usage     : $myReseqFrags->insert_gaps2frag($seqobj)
 Function  : inserts gaps to given fragment, by 
               1.) check if frag has deletion(s) and insert apropriate gaps
               2.) iterating over all positions in the max_ins_hash 
                   (respecting insertions, arisen by other frags)
 Returns   : locatable sequence object, with inserted gaps
 Args      : $seqobj : locatable sequence object


=cut

sub insert_gaps2frag() {

  my ($self, $seqobj) = @_;

  #$seqobj=insert_gaps2mtref($seq, $fragment, $fragment_span, $fragment_end, \%max_ins_hash, \%frags, $startpos, $endpos);
  
  my $myvalue = $seqobj->desc();
  my $startpos12=$self->_get_cont_no(\$myvalue);
  my $endpos12=$self->_get_cont_no(\$myvalue);
  $myvalue=$seqobj->display_id();
  my ($fragment, $fragment1, $fragment_no, $fragment_span, $fragment_end, $chip);
           
  if ($myvalue =~ /[:]/g) {
    $fragment=substr($myvalue,0,pos($myvalue)-1);
    $fragment1=$fragment;
    $chip=substr($myvalue,pos($myvalue));
    $fragment_no=$self->_get_cont_no(\$fragment1);
    $fragment_span=$self->_get_cont_no(\$fragment1);
    $fragment_end=$self->_get_cont_no(\$fragment1);
  } else {
    $fragment="unknown";
    $self->warn("Element with id ".$seqobj->id()." and description ".$seqobj->desc()." doesnt meet required format");
  }
  my $startpos = $fragment_span;
  my $endpos = $fragment_end;
  #my @frag_info_array=();
  #iterate over all positions of insertion of maximal length
  my %temp_hash;
  while ( my ($key, $value) = each %{$self->{_frags_hash}} ) {
    if ($key eq $fragment) {
      %temp_hash=%$value;
  #    while ( my ($subkey,$subvalue) = each %$value ) {
  #        @frag_info_array=@$subvalue;
  #        push(@frag_info_array,$subkey);
  #    }
    }
  }
  #insert gaps if there is/are deletion(s)
  my $offset=0;
  foreach my $mkey (sort{$a<=>$b} keys %temp_hash) {
    if ($temp_hash{$mkey}[0] eq "ins") {
      if (exists($self->{_max_ins_hash}{$mkey})) {
        $offset+=length($temp_hash{$mkey}[1]);
      }
    }
    if ($temp_hash{$mkey}[0] eq "del") {
      my $npos=$mkey-$startpos-12+$offset-1;
      Bio::SeqUtils->mutate($seqobj, Bio::LiveSeq::Mutation->new(
        -seq => "-",
        -pos => $npos,
        -len => 0));
    }
  }
  $offset=0;  
  my $gapsum=0;
  #fuer jede max_insertion testen ob sie im aktuellen fragment vorkommt, hash durchlaufen
  foreach my $maxkey (sort{$a<=>$b} keys %{$self->{_max_ins_hash}}) {
    my $flag=1;
    foreach my $mkey (keys %temp_hash) {
      if ($temp_hash{$mkey}[0] eq "ins") {
        #ist diese position ebenfalls im fragment definiert
        if ($mkey == $maxkey) {
          if (length($temp_hash{$mkey}[1]) < $self->{_max_ins_hash}{$maxkey} ) {
            my $gap=$self->_create_gap($self->{_max_ins_hash}{$maxkey}-length($temp_hash{$mkey}[1]));
            my $npos=$maxkey-$startpos+$gapsum+2-12+length($temp_hash{$mkey}[1]);
            Bio::SeqUtils->mutate($seqobj, Bio::LiveSeq::Mutation->new(
              -seq => $gap,
              -pos => $npos,
              -len => 0));
            $gapsum+=$self->{_max_ins_hash}{$maxkey};
          }
          elsif (length($temp_hash{$mkey}[1]) == $self->{_max_ins_hash}{$maxkey}) {
            $gapsum+=$self->{_max_ins_hash}{$maxkey};
          }
          $flag=0;          
        }
      }
    }
    #wenn position nicht im aktuellen fragment (temp_hash) definiert ist, dann schauen, ob die position innerhalb des fragments liegt
    #und ggf. eine luecke insertieren
    if ($flag) {
      #if ($)
      if ($maxkey >= $startpos+$startpos12-2 and $maxkey <= $startpos+$endpos12-2) {
      my $gap=$self->_create_gap($self->{_max_ins_hash}{$maxkey});
      my $npos=$maxkey-$startpos+$gapsum+2-12;
      foreach my $key (keys %{$self->{_refseq_max_ins_hash}}) {
        if ($fragment_span>$key) {$npos+=$self->{_refseq_max_ins_hash}{$key};}
      }      
      if ($npos<=$seqobj->length) {
      Bio::SeqUtils->mutate($seqobj, Bio::LiveSeq::Mutation->new(
        -seq => $gap,
        -pos => $npos,
        -len => 0));
      $gapsum+=$self->{_max_ins_hash}{$maxkey};
      }
      }
    }
  }
  my $new_id= sprintf '%.6d', $startpos;
  my $locseq = new Bio::LocatableSeq(
              -seq => $seqobj->seq,
              -id => $new_id."_".$seqobj->id,
              -start => $fragment_span-1, -end => $endpos12+$fragment_span-12);
  return $locseq;
}


=head2 insert_gaps2reference_sequence()

 Title     : insert_gaps2reference_sequence()
 Usage     : $myReseqFrags->insert_gaps2frag($seqobj)
 Function  : iterate over all positions in the pos_hash and insert gaps in the entire reference sequence
 Returns   : locatable sequence object, with inserted gaps
 Args      : $seqobj : locatable sequence object


=cut

sub insert_gaps2reference_sequence() {

  my ($self, $seqobj) = @_;
  my $offset=0;
  
  foreach my $key (keys %{$self->{_refseq_max_ins_hash}}) {
    $self->{_max_ins_hash}{$key}=$self->{_refseq_max_ins_hash}{$key};
  }
  foreach my $key (sort{$a<=>$b} keys %{$self->{_max_ins_hash}}) {
    my $gap=$self->_create_gap($self->{_max_ins_hash}{$key});
    Bio::SeqUtils->mutate($seqobj, Bio::LiveSeq::Mutation->new(
    -seq => $gap,
    -pos => $key+$offset-12+1, #+$mito_start_position_offset,
    -len => 0));
    $offset+=length($gap);
  }
  my $locseq = new Bio::LocatableSeq(
      -seq => $seqobj->seq,
      -id => "000000_".$seqobj->id(),
      -name => $seqobj->id(),
      -start => 0, -end => $seqobj->length());
  foreach my $key (keys %{$self->{_refseq_max_ins_hash}}) {
    #$self->{_max_ins_hash}{$key}=$self->{_refseq_max_ins_hash}{$key};
    delete($self->{_max_ins_hash}{$key});
  }
  return $locseq;
}



=head2 _create_gap()

 Title     : _create_gap()
 Usage     : private Function, dont call directly
 Function  : creates a gap (string of contiguous '-' chars) of given length
 Returns   : string
 Args      : $length : length of gap


=cut

sub _create_gap() {
  my ($self, $length) = @_;
  my $gapchar="-";
  my $i=0;
  my $gap="";
  while ($i<$length) {
    $i++;
    $gap.="-";
  }
  return $gap;
}




#$sub _print_H() {

#  my ($self,$hash_ref) = @_;
#  while ( my ($family, $roles) = each %$hash_ref ) {
#    print "$family: $roles\n";
#  }
#}




#sub _print_HoH() {
#
#  my ($self,$hash_ref) = @_;
#  #my %hash=%$hash_ref;
#  while ( my ($family, $roles) = each %$hash_ref ) {
#    print "$family: ";
#    while ( my ($role, $person) = each %{$roles} ) {
#      print "$role: $person\n";
#    }
#  }
#          
#}


#sub _print_HoHoA() {


#  my ($self,$member) = @_;
#  while ( my ($family, $roles) = each %{$self->{$member}} ) {
#    print "$family: ";
#    while ( my ($role, $person) = each %{$roles} ) {
#      print "$role: @$person\n";
#      #while ( my ($persons, $thing) = each %{$person} ) {
#      #  print "$persons=$thing";
#      #}
#    }
#  }
#}  

sub _print_base_hash() {

  my ($self, $hash_ref) = @_;
  my %hash=%$hash_ref;
  while ( my ($family, $roles) = each %hash ) {
    #print "$family: ";
    while ( my ($role, $person) = each %$roles ) {
      print "$person\n";
      #push(@stats_array, $person);
    }
  }
}


=head2 _check_oligo_positions()

 Title     : _check_oligo_positions()
 Usage     : private Function, dont call directly
 Function  : check if we are in a area which is covered by redundand fragments, if not, jump to the next such area
 Returns   : $position : position of the next cluster of fragments
 Args      : $pos : position in the sequence


=cut

sub _check_oligo_positions() {
  my ($self, $i) = @_;
  my $test=0;
  #foreach my $maxkey (sort{$a<=>$b} keys %{$self->{_max_ins_hash}}) {
  my $startpos_old="";
  my $endpos_old="";
  foreach my $startpos (sort{$a<=>$b} keys %{$self->{_oligos2calc_hash}}) {

    if ($endpos_old ne "") {
      if ($i>$endpos_old and $i<$startpos) {
        return $startpos;
      }
    }
    #if ($startpos < $i and $endpos > $i) {
    #  $test=1;
    #  last;
    #}
    $startpos_old=$startpos;
    $endpos_old=$self->{_oligos2calc_hash}{$startpos};
  }
  return $i;
}

=head2 _consider_context()

 Title     : _consider_context()
 Usage     : private Function, dont call directly
 Function  : check if there are n's up and down of the position
 Returns   : true if no ns are in the specified area else false
 Args      : $seq : locatable sequence object (redundand fragment)
             $pos : current position
             $options_hash : hash of options


=cut


sub _consider_context() {

  my ($self, $seq, $pos, $left_flank, $right_flank, $allowed_n) = @_;
  my $start_pos=1;
  my $end_pos=$seq->length()-1;
  my $count_left=0;
  my $count_right=0;
  my $count_gaps_right=0;
  my $count_gaps_left=0;

  if ($pos-$left_flank>1) {
    $start_pos=$pos-$left_flank;
    $count_gaps_left=($seq->subseq($start_pos, $pos-1) =~ tr/-//);
    if ($pos-$left_flank-$count_gaps_left>1) {
      $start_pos-=$count_gaps_left;
    }
    $count_left = ($seq->subseq($start_pos, $pos-1) =~ tr/n//);
    
  }
  if ($pos+$right_flank<$seq->length()-1) {
    $end_pos=$pos+$right_flank;
    $count_gaps_right=($seq->subseq($pos+1, $end_pos) =~ tr/-//);
    if ($pos+$right_flank+$count_gaps_right<$seq->length()-1) {
      $end_pos+=$count_gaps_right;
    }
    $count_right = ($seq->subseq($pos+1, $end_pos) =~ tr/n//)
  }

  if (($count_right+$count_left)>$allowed_n or $seq->subseq($pos, $pos) eq "n") {return 0}
  else {return 1}
}

=head2 _consider_context_wrapper()

 Title     : _consider_context_wrapper()
 Usage     : private Function, dont call directly
 Function  : wrapps call of _consider_context, and perform 2 calls  ('symmetrically' to both flanks (if length differ), if flanks have different length)
 Returns   : true if at least one of the call returns true else false
 Args      : $seq : locatable sequence object (redundand fragment)
             $pos : current position
             $options_hash : hash of options


=cut


sub _consider_context_wrapper() {

  my ($self, $seq, $pos, $options_hash, $ref_base) = @_;
  my $flank1;
  my $flank2;
  my $start_pos=1;
  my $end_pos=$seq->length()-1;
  #distinguish between parameter set between insertions and sub/del
  #if we have an insertion
  my $allowed_n;
  if ($ref_base eq '-') {
    $allowed_n=$options_hash->{allowed_n_in_flank_ins};
    $flank1=$options_hash->{flank_left_ins};
    $flank2=$options_hash->{flank_right_ins};
  }
  #else
  else {
    $allowed_n=$options_hash->{allowed_n_in_flank};
    $flank1=$options_hash->{flank_left};
    $flank2=$options_hash->{flank_right};
  }
  
  if ($flank1 != $flank2) {
    if ($self->_consider_context($seq, $pos, $flank1, $flank2, $allowed_n) or $self->_consider_context($seq, $pos, $flank2, $flank1, $allowed_n)) {return 1}
    else {return 0}
  }
  else {
    return $self->_consider_context($seq, $pos, $flank1, $flank2, $allowed_n);
  }
}


#=head2 _get_max_ins_length() (not used)
#
# Title     : _get_max_ins_length()
# Usage     : private Function, dont call directly
# Function  : iterate through hash of max. insertions and calcs the inseriton of maximal length
# Returns   : length of max insertion (integer)
# Args      : None
#
#
#=cut
#
#sub _get_max_ins_length() {
#  my $self=shift;
#  my $max_ins=0;
#  foreach my $key (sort{$a<=>$b} keys %{$self->{_max_ins_hash}}) {
#    if ($self->{_max_ins_hash}{$key} => $max_ins+1) {
#      $max_ins=$self->{_max_ins_hash}{$key};
#    }
#  }
#  return $max_ins;
#}

=head2 _swap_insertion()

 Title     : _swap_insertion()
 Usage     : private Function, dont call directly
 Function  : iterate over hash of max insertions and check if there is a deletion one position before 
             the current pos and the current position belongs to an insertion defined in max_ins_hash
 Returns   : true/false
 Args      : $pos : current position of the insertion


=cut


sub _swap_insertion() {

  my ($self, $pos, $sequence) = @_;
  my $swap=0;
  #print "\ninside swap_insertion ($pos)\n";
  #print "offset: ".$self->_get_start_pos($pos-12, '_frags_max_ins_hash', 1)."\n";
  my $pos1=$pos;
  my $offset=$self->_get_start_pos($pos-12);
  $pos+=12-$offset;
  #print "newpos: $pos\n";
  foreach my $key (sort{$a<=>$b} keys %{$self->{_max_ins_hash}}) {
    #print "\tstartgap: ".($key)."\tendgap: ";
    #print ($key+$self->{_max_ins_hash}{$key}-1);
    #print " vs ".($pos)."\n";
    #print "alskdjlaksjdlaksjdlaskdjlaskdj\n\n\n";
    if ($key+$self->{_max_ins_hash}{$key}-1>=$pos and $key+2<=$pos) {
      #print "can=swap==============\n";
      #print $key." vs ".($pos)."\t";
      #$swap=$pos-$key;
      my $isgap=1;
      my $i=1;
      while ($isgap) {
        #print "letzte 5 pos der biserigen seq: ".substr($sequence,-5)."\n";
        if (substr($sequence,0-$i) eq "-") {
          $swap++;
          $i++;
        }
        else {
          $isgap=0;
        }
      }
      #print "swap: (#positions): $swap\n";
      last;
    }
  }
  return $swap;
}

=head2 calc_sequence()

 Title     : calc_sequence()
 Usage     : $myReseqFrags->calc_sequence($aln, $options_hash [, $filename]);
 Function  : iterates over each position of the sequence (first element in the alignment) and change calls if sufficient evidence can be obtained by the alignment of redundant fragments
 Returns   : sequence (string)
 Args      : $als : locatable sequence object (redundand fragment)
             $options_hash : hash of options for redundant fragments analysis
             [$filename] : filename for outputting results

=cut

sub calc_sequence() {

  my ($self, $aln, $options_hash, $filename_rawrow) = @_;
  my $final_seq="";
  my $i=1;
  #$self->_print_H($self->{_max_ins_hash});
  #$self->_print_H($self->{_frags_max_ins_hash});
  #print "\nOpitions for redundant fragments analysis:\n";
  #$self->_print_H($options_hash);
  $self->{_inserted_bases_hash}=undef;;
  while ($i<=$aln->length()) {
    my $seq_no=1;
    #print "h";
    #if ($i%2000==1) {print "\n".$i."\n";}
    my $ref_base;
    my $help_base="x";
    my @base_array=();
    my $base;
    my $uni_votum=1;
    my $i_neu=$self->_check_oligo_positions($i);
    foreach my $seq ($aln->each_seq) {
      my $offset=$self->_get_start_pos($seq->start());
      if ($seq_no==1) {
        $ref_base=$seq->subseq($i_neu,$i_neu);
        if ($filename_rawrow) {
          open (RAWROW, ">>$filename_rawrow") or die "Cannot open file\n $!\n";
          print RAWROW "\n $i $ref_base: ";
          close(RAWROW);
        }
        #print "$i: $ref_base\n";
        if ($i_neu != $i) {
          $final_seq.=$seq->subseq($i,$i_neu-1);
        }
        $i=$i_neu;
      }
      else {
        if ($seq->start() < $i-$offset and $seq->end()+$offset+2+12 > $i) {
          #print "$i id: ".$seq->id()." ";#.$seq->seq()."\n";
          my $cleared_pos=$i-($seq->start()+$offset);
            if ($cleared_pos<=$seq->length()) {
              $base=$seq->subseq($cleared_pos,$cleared_pos);
              #print $seq->subseq($cleared_pos,$cleared_pos);
              if ($filename_rawrow) {
                open (RAWROW, ">>$filename_rawrow") or die "Cannot open file\n $!\n";
                print RAWROW "$base";
                close(RAWROW);
              }              
              #print $base;#."\n";
              
              my $help_var=$seq->id();
              if ($options_hash->{consider_context}==1) {
                if ($self->_consider_context_wrapper($seq, $cleared_pos, $options_hash, $ref_base)) {
                  push(@base_array, $base);
                }
              }
              else {
                push(@base_array, $base);
              }
            }
          #}
        
        }
        #remove no more needed sequences
        if (($i) > ($seq->end()+$offset)) {
          $aln->remove_seq($seq);
        }
        #finish iteration, if startpos of current sequence is outside of current position
        if ($i<$seq->start()+$offset) {
          #print "last: $seq_no\n";
          last;
        }
      }
      $seq_no++;
      #print "\n";
    }
    #if ($i>290 and $i<310) {print "$i: $ref_base\n";}
    if (@base_array>0) {
      if ($ref_base ne '-' and $options_hash->{include_main_sequence}) {
        push(@base_array, $base);
      }
      
      my $alignment_depth=@base_array;
      my $newbase="";
      my $vote=$self->_calc_stats(\@base_array, $options_hash, $ref_base);
      ##explore deletions
      my $arstr=join("",@base_array);
      my $count = ($arstr=~ tr/-//);
      
      if ($ref_base ne $vote) {
        my $swap=0;
        if ($vote ne "x" and $vote ne "n") {
        
          if ($ref_base eq "n" and $options_hash->{depth_n}-1<$alignment_depth) {
            #deletionen von N-calls
            if ($vote eq "-" and $options_hash->{deletions}==1 and ($options_hash->{depth_del}-1<$alignment_depth)) {
              $newbase=$vote;
            }
            #sonstige calls
            elsif ($vote ne "-") {
              $newbase=$vote;
            }
          }
          #insertionen
          elsif ($ref_base eq "-" and $options_hash->{insertions}==1 and $options_hash->{depth_ins}-1<$alignment_depth) {
            if ($options_hash->{swap_ins}) {
              $swap=$self->_swap_insertion($i, $final_seq);
            }
            $newbase=$vote;
          }
          #deletionen von gecallten basen
          elsif ($vote eq "-" and $options_hash->{deletions}==1 and $options_hash->{recall}==1 and ($options_hash->{depth_del}-1<$alignment_depth)){
            #if ($options_hash->{swap_ins}) {
            #  $swap=$self->_swap_insertion($i, $final_seq);
            #}
            $newbase=$vote;
          }
          #sonstige recalls
          elsif ($vote ne "-" and $options_hash->{recall}==1 and ($options_hash->{depth_recall}-1<$alignment_depth)) {
            $newbase=$vote;
          }
          #print @base_array;
          #print " (Tiefe: $alignment_depth) \t\t$i (";
          #print $i+10;
          #print ") $ref_base vs $vote => $newbase \n";
          if ($filename_rawrow) {
            open (RAWROW, ">>$filename_rawrow") or die "Cannot open file\n $!\n";
            print RAWROW " => ";
            print RAWROW @base_array;
            print RAWROW " (Depth: $alignment_depth) \t\t$i ";
            print RAWROW "$ref_base vs $vote => $newbase";
            close(RAWROW);
          }
        }
        if ($newbase ne "") {
          if ($swap>=1) {
              #print "\n\nswap: $swap\nvor: ".substr($final_seq,0,0-$swap)."\n\n";
              #print "mitte: ".$newbase."\n";
              #print "danach: ".substr($final_seq,0-$swap)."\n\n";
              $final_seq=substr($final_seq,0,0-$swap).$newbase.substr($final_seq,0-$swap);
              #break();
          }
          else {
              $final_seq.=$newbase;
          }
        } else {
          $final_seq.=$ref_base;
        }
      }
      else {
        $final_seq.=$ref_base;
      }
    }
    else {
      $final_seq.=$ref_base;
    }
    $i++;
  }
  return $final_seq; 
}



=head2 _calc_stats()

 Title     : _calc_stats()
 Usage     : private Function, dont call directly
 Function  : calcs statistics of letters in the given array by applying the Statistics::Frequency Package
 Returns   : base, which exceeds given threshhold otherwise 'x'
 Args      : $base_array : reference to array of bases (alignment slice)
             $options_hash : hash of options for redundant fragments analysis

=cut

sub _calc_stats() {

  my ($self, $array_ref, $options_hash, $ref_base) = @_;
  my @stats_array=@$array_ref;

  my $f1 = Statistics::Frequency->new;
  my $res="x";
  $f1->add_data( \@stats_array );
  $f1->remove_elements('n');
  my $threshold;
  if ($ref_base eq '-') {
    $threshold=$options_hash->{ins_threshold};
  }
  #else
  else {
    $threshold=$options_hash->{call_threshold};
  }
  #insertion
  my %freq = $f1->frequencies;
  my $max=0;
  my $mbase;
  if (!($f1->frequencies_max and $f1->frequencies_sum and $threshold)) {
    print "1 ".$f1->frequencies_max."\n";
    print "2 ".$f1->frequencies_sum."\n";
    print "3 ".$threshold."\n";
  }
  if ($f1->frequencies_sum>0) {
    for my $base (keys %freq) {
      if ($max<$freq{$base}) {
        $max=$freq{$base};
        $mbase=$base;
      }
    }
    if ($mbase eq "-") {
      $threshold=$options_hash->{del_threshold};
      #print "DELETION, na fein...\n";
    }
    #print " rel:".$f1->frequencies_max/$f1->frequencies_sum;
    #print "max: ".$f1->frequencies_max." ";
    #print "sum: ".$f1->frequencies_sum." ($threshold) \n";
    if ($f1->frequencies_max/$f1->frequencies_sum >= ($threshold/100)) {
      $res=$mbase;
    }
  }
  return $res;
}


=head2 write2fasta() 

 Title     : write2fasta()
 Usage     : $myReseqFrags->write2fasta($seq, $id, $filename, $gap)
 Function  : write the alignment of redundant fragments refering to entire sequence to xls file
 Args      : $seq : sequence as string to write to file
             $id : identifier of the sequence
             $filename : name of fasta file
             $gap : 1 gaps are written to file, 0 gaps are removed. if omitted, gaps are also removed
                                           
=cut


sub write2fasta() {

  my ($self, $string, $id, $filename, $gap) = @_;
  #print "write2fasta: $string, $id, $filename\n";
  if (!$gap or $gap == 0) {
    $string=~ s/-//g;
  }
  my $out = Bio::SeqIO->new(-file => ">>$filename" ,
                             -format => 'Fasta');
  my $locseq = new Bio::LocatableSeq(
      -seq => $string,
      -id => $id,
      -name => $id,
      -start => 0, -end => length($string));
  #$out->write_seq($locseq);
  open (FILE, ">>$filename");
  print FILE "> $id\n";
  print FILE "$string\n";
  close(FILE);
}

=head2 write_alignment2xls()

 Title     : write_alignment2xls()
 Usage     : $myReseqFrags->write_alignment2xls($aln, $myworkbook, $sheetname)
 Function  : write the alignment of redundant fragments refering to entire sequence to xls file
 Args      : $aln : reference to the alignment
             $workbook : reference of the xls workbook
             $id : name of the worksheet, where the alignment is written
             $startcol : number of first fragment, which is written to xls. default is 1

=cut

sub write_alignment2xls() {
  my ($self, $aln, $workbook, $ind_id, $startcol) = @_;
  
  my $worksheet = $workbook->add_worksheet($ind_id);
  
  my $j=0;
  my $seqno=0;
  my $start=1;
  my $i=2;
  if ($startcol) {
    $start=$startcol;
  }
  foreach my $seq ( $aln->each_seq() ) {
    #print $seq->id();
    $j++;
    
    if ($j>$start) {
    $worksheet->write(1, $j-$start, $seq->id());
    $i=2;
    }
    my $test_complete_seq=($seq->id =~ /human_mtDNA_RCRS/);
    if ($test_complete_seq==1) {}
    else {$test_complete_seq=0;}
    
    if ($j>$start) {
    #refseq schreiben
    if (!$test_complete_seq) {
      my $offset=$self->_get_start_pos($seq->start());
      #print "vals: ".$offset."======".$seq->start."============".$seq->start()."===\n";        
      while ($i<=$seq->length()+1) {
        my $zeile=$i+$seq->start()+$offset;
        my $wert=$seq->subseq($i-1,$i-1);
        $worksheet->write($zeile, $j-$start, $wert);
        $i+=1;
      }
    }
    }
    #fragmente schreiben
    else {
      my $offset=$self->_get_start_pos($seq->start());
      while ($i<=$seq->length()+1) {
        my $zeile=$i+$seq->start()+$offset;
        my $wert=$seq->subseq($i-1,$i-1);
        #$worksheet->write($zeile, $j-$start, $wert);
        $worksheet->write($zeile, 0, $wert);
        $i+=1;
      }
    }
    #$seqno++;
  }
}


1;