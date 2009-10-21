#------------------------------------------------------------------------------
# PACKAGE : Bio::Microarray::Tools::ReseqChip
# PURPOSE : Analyse additional probe oligonucleotides of Resequencing Chips
# AUTHOR  : Marian Thieme
# CREATED : 21.09.2007
# REVISION: 
# STATUS  : Beta
#
# For documentation, run this module through pod2html
#------------------------------------------------------------------------------


=head1 NAME

Bio::Microarray::Tools::ReseqChip - Class for analysing additional probe oligonucleotides of Resequencing Chips (for instance Affy Mitochip v2.0)

=head1 SYNOPSIS

 #see also ../../../t/Microarray/Tools/ReseqChip.t for the corresponding test script and a complete working example. This script uses parameter settings that 
 #that had been optimized with more than 100 mitochip samples.

 use Bio::Microarray::Tools::ReseqChip;
 
 # fasta file with reference sequence context useed for the Chip (MitoChip in that case)
 my $refseq_chip=...;
 #... and design file that describes the addional probes, grouped by consecutive probes covering 20-100 consecutive positions referred to the reference sequence
 my $Affy_frags_design_filename=...;
 #format of the design file
 $format='affy_mitochip_v2';
 # positions that are missing with respect to the reference sequence (rCRS - cambridge reference sequence) are going to be marked, so numbering with respect to the rCRS is conform
 my %ref_seq_max_ins_hash=(3106 => 1);
 
 my $reseqfragSample=Bio::Tools::ReseqChip->new(
    $Affy_frags_design_filename,
    $format,
    \%ref_seq_max_ins_hash,
    $refseq_chip);

 my $aln = new Bio::SimpleAlign();
 my $in  = Bio::SeqIO->new(-file => $Affy_reseq_sample_fasta_file,
                          -format => 'Fasta');

 my %options_hash=(
   include_main_sequence => 1,
   insertions => 1,
   deletions => 1,
   depth_ins => 1,
   depth_del => 9,
   depth => 1,
   consider_context => 1,
  flank_left => 10,
  flank_right => 10,
  allowed_n_in_flank => 0,
  flank_left_ins => 4,
  flank_right_ins => 4,
  allowed_n_in_flank_ins => 1,
  flank_size_weak => 1,
  call_threshold => 55,
  ins_threshold => 35,
  del_threshold => 75,
  swap_ins => 1);
                         
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
 my $new_sequence=$reseqfragSample->calc_sequence($aln, \%options_hash [,"output_file"]);


=head1 DESCRIPTION

This Software module aim to infer information of the addtional oligonucleotide probes, covering different known variants.
Oligonucleotide Array based Resequencing is done in the local context of a reference sequence. Every position in 
the genomic areas of interest is interrogated using 8 different 25-mer oligonucleotide probes (forward and reverse strand). 
Their middle base varies across the four possible bases, while the flanking regions are identical 
with the reference sequence or its reverse strand respectively. For genomic regions with known variability across individuals, 
additional probes were added to the chip. They interrogate postions in the neighborhood of polymorphisms not only in the local context
of the reference sequence but also in the context of its known variants.
This software (ReseqChip.pm) is tested to work with MitoChip v2.0 Data, manufactured by Affymetrix and the parser (MitoChipV2Parser)
reads the probe design file (Affy mtDNA_Design_Annotion.xls) wich describes the design of the probes.

The software approaches the problem in the following way:
1. An alignment of the addtional probes to the reference sequence is created (taking account for insertions/deletion)
2. Based on that alignment each position, which is covered by at least one additional probe is investigated to find a consensus call.

This is done indirectly by excluding those probes, which appear to be inadequate for the individual. An indication for 
inadaquacy is a local accumulation of N-calls. We investigate calls in neighborhoods of length K around
each sequence position in all available local context probes and count the number of N-calls in them. 
That menas, in addition to the call obtained using the references sequence base call we obtain data from all alternative 
local background probes that were available for the current position. All probes with more then maxN N-calls in the 
K-neighborhood are excluded. Because it may happen that different candidate bases occur we introduce to more parameters minP and minU.
If more then minP probes remain after filtering and more then minU percent of them call the base x,
were x is the most frequently called base, then x is included in the final sequence, otherwise the letter N is included.


Assumption:
Position of gaps (marked by the character "-") which are inserted in several fragments and in the reference sequence itself 
refer to the numbering given by the reference sequence, as defined in the design file. 
The reference sequence and the design file must be given as input parameters.
The parameter hash containing the different program options (%options_hash), specifying the explained parameter and some further 
options is provided by the user: we recommend to only specify the following hash-values that correnspod to the hash-keys (parameters):
 - flank_right (K-neighbourhood up and downstream of the genotype position in question)
 - call_threshold (minU)
 - depth (minP)
 - allowed_n_flank (maxN in K-neighbourhood)
for each substitutions, insertions and deletions. Note, for deletions only the parameters "depth" and "call_threshold" is applicable.

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

Marian Thieme (marian.thieme@gmail.com)

=head1 COPYRIGHT

Copyright (c) 2007-2009 Institute of Functional Genomics, University Regensburg, granted by Baygene. All Rights Reserved.

This module is free software; you can redistribute it and/or modify it under the same terms as Perl itself.





=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::Microarray::Tools::ReseqChip;
                                   
use strict;
use warnings;

use base qw(Bio::Root::Root);


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
             
             
 Returns   : Returns a new ReseqChip object
 
 Args      : $Affy_frags_design_filename (Affymetrix xls design file, 
             for instance: mtDNA_design_annotation_FINAL.xls for mitochondrial Genome)
             
             $format (only xls is available, because its the format which is delivered by Affymetrix)
             
             [$reseq_max_ins_hash, $refseq] (insertions as hash (pos1 => insertions length1, pos2 => insertions length2, ...) 
             for reference sequence and $refseq (Locatable Sequence Object))
            
 Membervars: frags_hash		- contains all variations described by the affy_design_annotation file 
                                  (fragment_id => (pos1 => [muttype, mut, start, stop]), ... )
             max_ins_hash	- contains all (maximal) insertion, which are needed to build "alignable" 
                                  fragments by inserting appropriate gaps
             oligos2calc_hash	- describes the coverage of the redundand fragments (start => stop)
             refseq		- Locatable Sequence holds reference sequence
             refseq_max_ins_hash- insertions need to be applied to the reference sequence
             gapped_refseq	- reference sequence when insertions are done


=cut


sub new {

  #my ($class, $design_file_name, $format, $refseq_max_ins_hash, $refseq) = @_;
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($design_file_name, $format, $refseq_max_ins_hash, $refseq)=$self->_rearrange([qw(AFFY_DESIGN_FILENAME FORMAT_OF_DESIGN_FILE MAX_INSERTION_HASH_REFERENCE REFERENCE_SEQUENCE)], @args);

  $self->{_frags_hash}=undef;
  $self->{_max_ins_hash}=undef;
  $self->{_refseq}=undef;
  $self->{_gapped_refseq}=undef;
  $self->{_oligos2calc_hash}=undef;
  $self->{_refseq_max_ins_hash}=undef;
  $self->{_frags_max_ins_hash}=undef;

  ##
  $self->{_oligo_flank_length}=12;
  $self->throw("Must provide filename as first argument !") unless $design_file_name;

  $self->throw("Must specify format (only xls at present) as second argument !") unless $format;
  
  $self->warn("No reference sequence given !") unless $refseq;
  
  my %max_ins_hash=();
  
  
  if ($format eq "affy_mitochip_v2") {
    my $parser=Bio::Microarray::Tools::MitoChipV2Parser->new($design_file_name);
    $self->{_frags_hash}=$parser->{_frags_hash};
    $self->{_oligos2calc_hash}=$parser->{_oligos2calc_hash};
    $self->debug("Created array of redundant fragments. \n");
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

sub _print_frags_hash() {
	
	my ($self) = @_;
	my $correction=0;
	foreach my $key (sort{$a cmp $b} keys %{$self->{_frags_hash}}) {
		#print "$key: ";
		foreach my $subkey (sort{$a<=>$b} keys %{$self->{_frags_hash}{$key}}) {
    		my $subvalue=$self->{_frags_hash}{$key}{$subkey};
      		my @array= @{$self->{_frags_hash}->{$key}{$subkey}};
		    $self->debug("$key\t$subkey\t@array \n");
    		if (!(exists($self->{_max_ins_hash}{$subkey+$correction}))) {
        		if (@$subvalue[0] eq "ins") {
          			$self->{_max_ins_hash}{$subkey+$correction}=length(@$subvalue[1]);
        		}
      		} elsif ( $self->{_max_ins_hash}{$subkey+$correction} < length(@$subvalue[1]) ) {
  	        	$self->{_max_ins_hash}{$subkey+$correction}=length(@$subvalue[1]);
    	  	}
    	}
  	}	
	
}


=head2 _get_cont_no()

 Title     : _get_cont_no()
 Usage     : private Function, dont call directly
 Function  : finds first contiguous number in given string and cut that part away from the string
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
 Returns   : cummulative offset
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
    #print "inside get_start_pos\n";
  }
  my $offset=0; 
  foreach my $key (sort{$a<=>$b} keys %{$self->{$hashname}}) {
    if ($test) {
      #print "$key: ".$self->{$hashname}{$key}."\n";
    }
    if ($pos>$key-$self->{_oligo_flank_length}) {
    #if ($pos>=$key) {
      $offset+=$self->{$hashname}{$key};
      if ($test) {
        #print "($offset)\n";
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
      my $npos=$mkey-$startpos-$self->{_oligo_flank_length}+$offset+1;
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
  	my $roffset=0;
  	foreach my $rkey (sort{$a<=>$b} keys %{$self->{_refseq_max_ins_hash}}) {
  		if ($rkey<$maxkey) {
  			$roffset+=$self->{_refseq_max_ins_hash}{$rkey}
  		}
  	}
    my $flag=1;
    foreach my $mkey (keys %temp_hash) {
      if ($temp_hash{$mkey}[0] eq "ins") {
        #ist diese position ebenfalls im fragment definiert
        if ($mkey == $maxkey) {
          if (length($temp_hash{$mkey}[1]) < $self->{_max_ins_hash}{$maxkey} ) {
            my $gap=$self->_create_gap($self->{_max_ins_hash}{$maxkey}-length($temp_hash{$mkey}[1]));
            my $npos=$maxkey-$startpos+$gapsum+2-$self->{_oligo_flank_length}+length($temp_hash{$mkey}[1])-$roffset;
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
      my $npos=$maxkey-$startpos+$gapsum+2-$self->{_oligo_flank_length}-$roffset;
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
              -start => $fragment_span-1, -end => $endpos12+$fragment_span-$self->{_oligo_flank_length}-1);
  return $locseq;
}


=head2 insert_gaps2reference_sequence()

 Title     : insert_gaps2reference_sequence()
 Usage     : $myReseqFrags->insert_gaps2frag($seqobj)
 Function  : iterate over all positions in the pos_hash and insert gaps in the reference sequence
 Returns   : locatable sequence object, with inserted gaps
 Args      : $seqobj : locatable sequence object


=cut

sub insert_gaps2reference_sequence() {

  my ($self, $seqobj) = @_;
  my $offset=0;
  
  foreach my $key (keys %{$self->{_refseq_max_ins_hash}}) {
    #$self->{_max_ins_hash}{$key}=$self->{_refseq_max_ins_hash}{$key};
    my $gap=$self->_create_gap($self->{_refseq_max_ins_hash}{$key});
    Bio::SeqUtils->mutate($seqobj, Bio::LiveSeq::Mutation->new(
    -seq => $gap,
    -pos => $key-$self->{_oligo_flank_length}+1, #+$mito_start_position_offset,
    -len => 0));
    #$offset+=length($gap);
  }
  foreach my $key (sort{$a<=>$b} keys %{$self->{_max_ins_hash}}) {
    my $gap=$self->_create_gap($self->{_max_ins_hash}{$key});
    Bio::SeqUtils->mutate($seqobj, Bio::LiveSeq::Mutation->new(
    -seq => $gap,
    -pos => $key+$offset-$self->{_oligo_flank_length}+1, #+$mito_start_position_offset,
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
 Returns   : string, representing the gaps
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



sub _print_base_hash() {

  my ($self, $hash_ref) = @_;
  my %hash=%$hash_ref;
  while ( my ($family, $roles) = each %hash ) {
    #print "$family: ";
    while ( my ($role, $person) = each %$roles ) {
      #print "$person\n";
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
  my $startpos_old="";
  my $endpos_old="";
  my $lastone=1;
  foreach my $startpos (sort{$a<=>$b} keys %{$self->{_oligos2calc_hash}}) {
    my $endpos=$self->{_oligos2calc_hash}{$startpos};  
    #if ($i>16239) {print("$i check oligo region: $startpos $endpos if:".($i>$endpos_old and $i<$startpos)."else:".($i>=$startpos and $i<=$endpos)."\n")}
    if ($endpos_old ne "") {
      if ($i>$endpos_old and $i<$startpos) {
        return $startpos;
      } elsif ($i>=$startpos and $i<=$endpos) {
        $lastone=1;
        last;
      }
    }
    ##else i might be before the first cluster
    else {
      if ($startpos > $i) {
        return $startpos;
      }
    }

    $startpos_old=$startpos;
    $endpos_old=$endpos;
    if ($i>$endpos_old) {
      $lastone=0
    }
  }
  if ($lastone==0) {
    return 0;
  }
  return $i;
}

=head2 _consider_context()

 Title     : _consider_context()
 Usage     : private Function, dont call directly
 Function  : check if there are n's up and down of the position
 Returns   : true if no ns are in the specified area else false
 Args      : $seq : locatable sequence object (reference sequenc or additional probe)
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
    if ($start_pos<$pos) {
      $count_gaps_left=($seq->subseq($start_pos, $pos-1) =~ tr/-//);
      if ($pos-$left_flank-$count_gaps_left>1) {
        $start_pos-=$count_gaps_left;
      }
      $count_left = ($seq->subseq($start_pos, $pos-1) =~ tr/n//);
    }
  }
  if ($pos+$right_flank<$seq->length()-1) {
    $end_pos=$pos+$right_flank;
    if ($pos<$end_pos) {
      $count_gaps_right=($seq->subseq($pos+1, $end_pos) =~ tr/-//);
      if ($pos+$right_flank+$count_gaps_right<$seq->length()-1) {
        $end_pos+=$count_gaps_right;
      }
      $count_right = ($seq->subseq($pos+1, $end_pos) =~ tr/n//)
    }
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
    #check if ends are within the flank size
    if ($options_hash->{flank_size_weak}) {
      if ($pos+$flank2>$seq->length()) {
        $flank2=$seq->length()-$pos;
      }
      if ($pos-$flank1<1) {
        $flank1=$pos;
      }
    }
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
  my $pos1=$pos;
  my $offset=$self->_get_start_pos($pos-$self->{_oligo_flank_length});
  $pos+=$self->{_oligo_flank_length}-$offset;
  foreach my $key (sort{$a<=>$b} keys %{$self->{_max_ins_hash}}) {
    if ($key+$self->{_max_ins_hash}{$key}-1>=$pos and $key+2<=$pos) {
      my $isgap=1;
      my $i=1;
      while ($isgap) {
        if (substr($sequence,0-$i) eq "-") {
          $swap++;
          $i++;
        }
        else {
          $isgap=0;
        }
      }
      last;
    }
  }
  return $swap;
}

=head2 calc_sequence()

 Title     : calc_sequence()
 Usage     : $myReseqFrags->calc_sequence($aln, $options_hash [, $filename]);
 Function  : iterates over each position of the sequence (first element in the alignment) and changes calls if sufficient evidence can be obtained by the alignment of redundant fragments
 Returns   : sequence (string)
 Args      : $als : locatable sequence object (redundand fragment)
             $options_hash : hash of options for redundant fragments analysis
             [$filename] : filename for outputting results

=cut

sub calc_sequence() {

  my ($self, $aln, $options_hash, $filename_rawrow) = @_;

  my $final_seq="";
  my $start_c=1;
  my $stop_c=0;
  my $chip_offset=0;
  my $range=0;
  if ($options_hash->{start_pos}) {
    $start_c=$options_hash->{start_pos}
  }
  if ($options_hash->{stop_pos}) {
    $stop_c=$options_hash->{stop_pos}
  }
  my $i=$start_c;
  
  my $stop=0;
  my $output_rawrow="";
  my $no_removed=0;
  my $index_ex=0;
  if ($filename_rawrow) {
    my $seq = $aln->get_seq_by_pos(1);
    open (my $RAWROW, ">>$filename_rawrow") or die "Cannot open file\n $!\n";
    print $RAWROW "\n".$seq->id;
    close($RAWROW);
  }
  while ($i<=$aln->length()) {
    
    my $seq_no=1;
    my $ref_base;
    my $help_base="x";
    my @base_array=();
    my $uni_votum=1;
    ##check if there are fragments for current position i
    my $i_neu=$self->_check_oligo_positions($i);
    if ($i_neu==0) {
      $i_neu=$aln->length();
    }

    my $count=0;
    my $output_rawrow_tmp="";
    my $not_only_ref=0;
    my $alt_seq;
        
    foreach my $seq ($aln->each_seq) {
      if ($options_hash->{alternative_sequence_hash}) {
        $alt_seq=$options_hash->{alternative_sequence_hash}{$seq->id()};
      }
      my $offset=$self->_get_start_pos($seq->start());
      if ($seq_no==1) {
        $ref_base=$seq->subseq($i_neu,$i_neu);
        if ($filename_rawrow) {
          $output_rawrow_tmp.= "\n".($i_neu+$chip_offset)." ($ref_base) : ";
        }
        if ($i>$stop_c and $stop_c>0) {
          $stop=1;
          $final_seq.=$seq->subseq($i,$seq->length());
          last;
        }
        if ($i_neu != $i) {
          $final_seq.=$seq->subseq($i,$i_neu-1);
        }
        $i=$i_neu;
      }
      
      ##add base to basearray if it fullfill the criteria
      ##differ between alternative base calls (in case of an insertion)
      if ($ref_base eq "-" and $options_hash->{alternative_sequence_hash}) {
        ($not_only_ref, $count, $output_rawrow_tmp)=$self->_augment_base_array($alt_seq, $ref_base, \@base_array, $not_only_ref, $count, $i, $offset, $seq_no, $options_hash, $filename_rawrow, $output_rawrow_tmp);
      ##and "normal" base calls for the alternative probes
      } else {
        ($not_only_ref, $count, $output_rawrow_tmp)=$self->_augment_base_array($seq, $ref_base, \@base_array, $not_only_ref, $count, $i, $offset, $seq_no, $options_hash, $filename_rawrow, $output_rawrow_tmp);
      }
      
      #remove no more needed sequences
      if (($i) >= ($seq->end()+$offset+2+$self->{_oligo_flank_length})) {
        $aln->remove_seq($seq);
      }
      #finish iteration, if startpos of current fragment/sequence is outside of current position in the alignment
      if ($i<$seq->start()+$offset) {
        last;
      }
      $seq_no++;
    }
    
    if ($stop) {
      last;
    }
    #at least one nonref base is available
    if ($not_only_ref) {
      $output_rawrow.=$output_rawrow_tmp;
      ($final_seq, $output_rawrow)=$self->_get_consensus_call($ref_base, \@base_array, $count, $final_seq, $options_hash, $filename_rawrow, $output_rawrow, $i);
    } else {
	      $final_seq.=$ref_base;
    }

    $i++;
  }
  if ($filename_rawrow) {
    open (my $RAWROW, ">>$filename_rawrow") or die "Cannot open file\n $!\n";
    print $RAWROW $output_rawrow;
    close($RAWROW);
  }

  return $final_seq; 
}



=head2 _augment_base_array()

 Title     : _augment_base_array()
 Usage     : dont call directly
 Function  : pushs base into the array of candidate bases, if the base meets criteria (options_hash->{allowed_n_in_flank}, options_hash->{flank_left/right})
 Returns   : 3 values:
             $not_only_ref 0/1: whether the array of candidate bases contains only reference base
             $count integer: number of candidate base currently in the array
             $output_rawrow_tmp string: log message, currently not used.
 Args      : lot of args, see the function itself

=cut

sub _augment_base_array() {

      my ($self, $seq, $ref_base, $base_array, $not_only_ref, $count, $i, $offset, $seq_no, $options_hash, $filename_rawrow, $output_rawrow_tmp) = @_;
      if ($seq->start() < $i-$offset and ($seq->end()+$offset+$self->{_oligo_flank_length}) > $i) {
        my $cleared_pos=$i-($seq->start()+$offset);
        if ($cleared_pos<=$seq->length()) {
            my $base=$seq->subseq($cleared_pos, $cleared_pos);

            if ($base ne 'n') {
              if ($filename_rawrow) {
                $output_rawrow_tmp.= $base;
                $count++;
              }              
              
              my $help_var=$seq->id();
              if ($options_hash->{consider_context}==1) {
                if ($self->_consider_context_wrapper($seq, $cleared_pos, $options_hash, $ref_base)) {
                  if ( $seq_no==1) {
                    if ($options_hash->{include_main_sequence}) {
                      push(@$base_array, $base);
                    }
                  } else {
                    push(@$base_array, $base);
                    $not_only_ref=1;
                  }
                }
              } else {
                push(@$base_array, $base);
              }
            }
        }
      }
      return($not_only_ref, $count, $output_rawrow_tmp);

}


=head2 _get_consensus_call()

 Title     : _get_consensus_call()
 Usage     : dont call directly
 Function  : based on the array of candidate base a consensus base is determined, possibly the n-call letter n is introduced
 Returns   : 2 values:
             $final_seq string: sequence processed so far
             $output_rawrow_tmp string: log message.
 Args      : lot of args, see the function itself

=cut

sub _get_consensus_call() {

      my ($self, $ref_base, $base_array, $count, $final_seq, $options_hash, $filename_rawrow, $output_rawrow, $pos) = @_;
    
      my $alignment_depth=@$base_array;
      my $newbase="";
      my $vote=$self->_calc_stats($base_array, $options_hash, $ref_base, $pos);
      my $arstr=join("", @$base_array);
      if ($ref_base ne $vote) {
        my $swap=0;
        if ($filename_rawrow) {
          $output_rawrow.= " ($count) -> $arstr ($alignment_depth)";
        }
        if ($vote ne "x" and $vote ne "n") {
          ###deletions
          if ($options_hash->{deletions}==1) {
            if ( $vote eq "-" and $ref_base ne "-" and $options_hash->{depth_del}<=$alignment_depth ) {
              $newbase=$vote;
            }
          }
          ###insertions
          if ($options_hash->{insertions}==1) {
             if ($vote ne "-" and $ref_base eq "-" and $options_hash->{depth_ins}<=$alignment_depth ) {
               $newbase=$vote;
            }
          }
          ###substitutions
          if ($vote ne "-" and $ref_base ne "-") {
            if ($options_hash->{depth}<=$alignment_depth) {
              $newbase=$vote;
            }
          }
          if ($filename_rawrow) {
            $output_rawrow.= "\t$ref_base vs $vote => $newbase";
          }
        } else {
          if ($options_hash->{call_n}) {
            $newbase=$vote;
          } else {
            $newbase=$ref_base;
          }
          if ($filename_rawrow) {
            $output_rawrow.= "\t$ref_base vs $vote => $newbase (no consensus call)";
          }

        }
        
        if ($newbase ne "") {
          if ($swap>=1) {
              $final_seq=substr($final_seq,0,0-$swap).$newbase.substr($final_seq,0-$swap);
          }
          else {
              $final_seq.=$newbase;
          }
        } else {
          if ($options_hash->{call_n}) {
            $final_seq.="n";
          } else {
            $final_seq.=$ref_base;
          }
        }
      }
      else {
        $final_seq.=$ref_base;
        if ($filename_rawrow) {
          $output_rawrow.= "\tnothing to do (vote == ref)";
        }
      }
      return($final_seq, $output_rawrow);
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

  my ($self, $array_ref, $options_hash, $ref_base, $pos) = @_;
  my @stats_array=@$array_ref;

  my $f1 = Statistics::Frequency->new;
  my $res="x";
  $f1->add_data( \@stats_array );
  $f1->remove_elements('n');
  my $threshold;
  if ($ref_base eq '-') {
    $threshold=$options_hash->{ins_threshold};
  } else {
    $threshold=$options_hash->{call_threshold};
  }
  my %freq = $f1->frequencies;
  my $max=0;
  my $mbase;

  ##get most frequent base
  if ($f1->frequencies_sum>0) {
    for my $base (keys %freq) {
      if ($max<$freq{$base}) {
        $max=$freq{$base};
        $mbase=$base;
      }
    }
    if ($mbase eq "-" and $ref_base ne "-") {
      $threshold=$options_hash->{del_threshold};
    }
    if ($f1->frequencies_max/$f1->frequencies_sum >= ($threshold/100)) {
      $res=$mbase;
    } elsif ($options_hash->{call_n}) {
    	$res="n";
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
  $out->write_seq($locseq);
  #open (my $FILE, ">>$filename");
  #print $FILE "> $id\n";
  #print $FILE "$string\n";
  #close($FILE);
}

=head2 write_alignment2xls()

 Title     : write_alignment2xls()
 Usage     : $myReseqFrags->write_alignment2xls($aln, $myworkbook, $sheetname)
 Function  : write the alignment ($aln object) to xls file.
 Args      : $aln : reference to the alignment
             $workbook : reference of the xls workbook
             $id : name of the worksheet, where the alignment is written
             $startcol : number of first fragment, which is written to xls. default is 1
             $refid : identifier for the reference sequence, so that this sequence can written to the leftmost column

=cut

sub write_alignment2xls() {
  #my ($self, $aln, $workbook, $ind_id, $startcol, $refid) = @_;
  my ($self, $aln, $workbook, $ind_id, $refid, $startcol) = @_;
  #print "Refid: $refid\n";
  #my $sheetname=$ind_id=~s/human_mtDNA_RCRS://g;
  my $sheetname=$ind_id=~s/$refid://g;
  my $worksheet = $workbook->add_worksheet($ind_id);
  
  my $j=0;
  my $seqno=0;
  my $start=1;
  my $i=2;
  my $k=1;
  if ($startcol) {
    $start=$startcol;
  }
  my $totlength=0;
  my $totlength_crs=0;
  my $totnumber=0;
  foreach my $seq ( $aln->each_seq() ) {
    #print $seq->id();
    $j++;
    $totnumber++;
    if ($j>$start) {
      $worksheet->write(1, $j-$start+1, $seq->id());
      $i=2;
    }
    my $test_complete_seq=($seq->id =~ /$refid/);
    #my $test_complete_seq=($seq->id =~ /human_mtDNA_RCRS/);
    if ($test_complete_seq) {$test_complete_seq=1}
    else {$test_complete_seq=0;}
    
    if ($j>$start) {
      #alternative fragments
      if (!$test_complete_seq) {
        $totlength+=$seq->length();
        my $offset=$self->_get_start_pos($seq->start());
        #print "vals: ".$offset."======".$seq->start."============".$seq->start()."===\n";
        while ($i<=$seq->length()+1) {
          my $zeile=$i+$seq->start()+$offset;
          my $wert=$seq->subseq($i-1,$i-1);
          $worksheet->write($zeile, $j-$start+1, $wert);
          $i+=1;
        }
      }
    }
    #refseq
    else {
      #print "Ref: ".$seq->length()."\n";
      $totlength_crs+=$seq->length();
      my $offset=$self->_get_start_pos($seq->start());
      while ($i<=$seq->length()+1) {
        my $zeile=$i+$seq->start()+$offset;
        my $wert=$seq->subseq($i-1,$i-1);
        #$worksheet->write($zeile, $j-$start, $wert);
        if ($wert ne "-" or exists($self->{_refseq_max_ins_hash}{$k})) {
        	
	        $worksheet->write($zeile, 0, $k+$self->{_oligo_flank_length});
	        $k++;
        }
        #if (exists($self->{_refseq_max_ins_hash}{$k}) ) {
       # 	print("exists: ($k)".$self->{_refseq_max_ins_hash}{$k});
        #	$k+=$self->{_refseq_max_ins_hash}{$k};
        #	print("und jetzt: $k\n");
        #	
        #}
        
        
        $worksheet->write($zeile, 1, $wert);
        $i+=1;
      }
    }
    #$seqno++;
  }
  #print "$totlength $totnumber $totlength_crs\n";
}


1;