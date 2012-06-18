#
# bioperl module for Bio::LiveSeq::Transcript
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::LiveSeq::Transcript - Transcript class for LiveSeq

=head1 SYNOPSIS

  # documentation needed

=head1 DESCRIPTION

This stores information about coding sequences (CDS).
The implementation is that a Transcript object accesses a collection of
Exon objects, inferring from them the nucleotide structure and sequence.

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::Transcript;

use strict;
# use Carp qw(carp cluck);
use Bio::LiveSeq::Exon; # uses Exon to create new exon in case of deletion
use base qw(Bio::LiveSeq::SeqI);

=head2 new

  Title   : new
  Usage   : $transcript = Bio::LiveSeq::Transcript->new(-exons => \@obj_refs);

  Function: generates a new Bio::LiveSeq::Transcript
  Returns : reference to a new object of class Transcript
  Errorcode -1
  Args    : reference to an array of Exon object references

=cut

sub new {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my ($obj,%transcript);

  my @exons=@{$args{-exons}};

  $obj = \%transcript;
  $obj = bless $obj, $class;

  unless (@exons) {
    $obj->warn("$class not initialised because exons array empty");
    return(-1);
  }

  # now useless, after start and end methods have been overridden here
  my $firstexon = $exons[0];
  #my $lastexon = $exons[-1];
  #my $start = $firstexon->start;
  #my $end = $lastexon->end;
  my $strand = $firstexon->strand;
  my $seq = $firstexon->{'seq'};
  $obj->alphabet('rna');

  unless (_checkexons(\@exons)) {
    $obj->warn("$class not initialised because of problems in the exon structure");
    return(-1);
  }
  $obj->{'strand'}=$strand;
  $obj->{'exons'}=\@exons;
  $obj->{'seq'}=$seq;

  # set Transcript into each Exon
  my $exon;
  foreach $exon (@exons) {
    $exon->{'transcript'}=$obj;
  }
  return $obj;
}


=head2 all_Exons

 Title   : all_Exons
 Usage   : $transcript_obj->all_Exons()
 Function: returns references to all Exon objects the Transcript is composed of
 Example : foreach $exon ($transcript->all_Exons()) { do_something }
 Returns : array of object references
 Args    : none

=cut

sub all_Exons {
  my $self=shift;
  my $exonsref=$self->{'exons'};
  my @exons=@{$exonsref};
  my @newexons;
  my $exon;
  foreach $exon (@exons) {
    unless ($exon->obj_valid) {
      $self->warn("$exon no more valid, start or end label lost, skipping....",1); # ignorable
    } else {
      push(@newexons,$exon);
    }
  }
  if ($#exons != $#newexons) {
    # update exons field
    $self->{'exons'}=\@newexons;
  }
  return (@newexons);
}

=head2 downstream_seq

 Title   : downstream_seq
 Usage   : $transcript_obj->downstream_seq()
         : $transcript_obj->downstream_seq(64)
 Function: returns a string of nucleotides downstream of the end of the
           CDS. If there is some information of the real mRNA, from features in
           an attached Gene object, it will return up to those boundaries.
           Otherwise it will return 1000 nucleotides.
           If an argument is given it will override the default 1000 number
           and return instead /that/ requested number of nucleotides.
           But if a Gene object is attached, this argument will be ignored.
 Returns : string
 Args    : an optional integer number of nucleotides to be returned instead of
           the default if no gene attached

=cut

sub downstream_seq {
  my ($self,$howmany)=@_;
  my $str;
  if (defined ($howmany)) {
    unless ($howmany > 0) {
      $self->throw("No sense in asking less than 1 downstream nucleotides!");
    }
  } else {
    unless ($self->{'seq'}->alphabet eq 'rna') { # if rna retrieve until the end
      #$str=$DNAobj->labelsubseq($self->end,undef,undef,"unsecuremoderequested");
      #return(substr($str,1)); # delete first nucleotide that is the last of Transcript
      if ($self->gene) { # if there is Gene object attached fetch relevant info
	$str=$self->{'seq'}->labelsubseq($self->end,undef,$self->gene->maxtranscript->end); # retrieve from end of this Transcript to end of the maxtranscript
	$str=substr($str,1); # delete first nucleotide that is the last of Transcript
	if (CORE::length($str) > 0) {
	  return($str);
	} else { # if there was no downstream through the gene's maxtranscript, go the usual way
	  $howmany = 1000;
	}
      } else {
	$howmany = 1000;
      }
    }
  }
  my @exons=$self->all_Exons;
  my $strand=$self->strand();
  my $lastexon=$exons[-1];
  my $lastexonlength=$lastexon->length;
  # $howmany nucs after end of last exon
  #my $downstream_seq=$lastexon->subseq($lastexonlength+1,undef,$howmany);
  my $downstream_seq;

  if ($howmany) {
      $downstream_seq=substr($lastexon->labelsubseq($self->end,$howmany,undef,"unsecuremoderequested"),1);
  } else {
    if ($strand == 1) {
      $downstream_seq=substr($lastexon->labelsubseq($self->end,undef,$self->{'seq'}->end,"unsecuremoderequested"),1);
    } else {
      $downstream_seq=substr($lastexon->labelsubseq($self->end,undef,$self->{'seq'}->start,"unsecuremoderequested"),1);
    }
  }
  return $downstream_seq;
}

=head2 upstream_seq

 Title   : upstream_seq
 Usage   : $transcript_obj->upstream_seq()
         : $transcript_obj->upstream_seq(64)
 Function: just like downstream_seq but returns nucleotides before the ATG
 Note    : the default, if no Gene information present and no nucleotides
           number given, is to return up to 400 nucleotides.

=cut

sub upstream_seq {
  my ($self,$howmany)=@_;
  if (defined ($howmany)) {
    unless ($howmany > 0) {
      $self->throw("No sense in asking less than 1 upstream nucleotides!");
    }
  } else {
    unless ($self->{'seq'}->alphabet eq 'rna') { # if rna retrieve from the start
      if ($self->gene) { # if there is Gene object attached fetch relevant info
	my $str=$self->{'seq'}->labelsubseq($self->gene->maxtranscript->start,undef,$self->start); # retrieve from start of maxtranscript to start of this Transcript
	chop $str; # delete last nucleotide that is the A of starting ATG
	if (length($str) > 0) {
	  return($str);
	} else { # if there was no upstream through the gene's maxtranscript, go the usual way
	  $howmany = 400;
	}
      } else {
	$howmany = 400;
      }
    }
  }
  my @exons=$self->all_Exons;
  my $firstexon=$exons[0];
  
  my $upstream_seq;
  my $strand=$self->strand();

  if ($howmany) {# $howmany nucs before begin of first exon
    my $labelbefore=$firstexon->label(-$howmany,$firstexon->start);
    if ($labelbefore < 1) {
      if ($strand == 1) {
	$labelbefore=$self->{'seq'}->start;
      } else {
	$labelbefore=$self->{'seq'}->end;
      }
    }
    $upstream_seq=$firstexon->labelsubseq($labelbefore,undef,$firstexon->start,"unsecuremoderequested");
    chop $upstream_seq;
  } else {
    if ($strand == 1) {
      $upstream_seq=$firstexon->labelsubseq($self->{'seq'}->start,undef,$self->start,"unsecuremoderequested");
      chop $upstream_seq; # delete last nucleotide that is the A of starting ATG
    } else {
      $upstream_seq=$firstexon->labelsubseq($self->{'seq'}->end,undef,$self->start,"unsecuremoderequested");
      chop $upstream_seq; # delete last nucleotide that is the A of starting ATG
    }
  }
  return $upstream_seq;
}

# These get redefined here, overriding the SeqI one because they draw their
# information from the Exons a Transcript is built of
# optional argument: firstlabel. If not given, it checks coordinate_start
#                                This is useful when called by Translation
#                                also used by _delete
sub label {
  my ($self,$position,$firstlabel)=@_;
  unless ($position) {  # if position = 0 complain ?
    $self->warn("Position not given or position 0");
    return (-1);
  }
  my ($start,$end,$strand)=($self->start(),$self->end(),$self->strand());
  my ($label,@labels,$length,$arraypos);
  unless (defined ($firstlabel)) {
    $firstlabel=$self->coordinate_start; # this is inside Transcript obj
  }
  my $coord_pos=$self->_inside_position($firstlabel);
  $length=$self->length;
  #if ($strand == 1) {
    if ($position < 1) {
      $position++; # to account for missing of 0 position
    }
    $arraypos=$position+$coord_pos-2;
    #print "\n=-=-=-=-DEBUG: arraypos $arraypos, pos $position, coordpos: $coord_pos";
    if ($arraypos < 0) {
      $label=$self->{'seq'}->label($arraypos,$start,$strand); #?
    } elsif ($arraypos >= $length) {
      $label=$self->{'seq'}->label($arraypos-$length+2,$end,$strand); #?
    } else { # inside the Transcript
      @labels=$self->all_labels;
      $label=$labels[$arraypos];
    }
  #}
}

# argument: label
# returns: position of label according to coord_start
# errorcode: 0 label not found
# optional argument: firstlabel. If not given, it checks coordinate_start
#                                This is useful when called by Translation
sub position {
  my ($self,$label,$firstlabel)=@_;
  unless ($self->{'seq'}->valid($label)) {
    $self->warn("label is not valid");
    return (0);
  }
  unless (defined ($firstlabel)) {
    $firstlabel=$self->coordinate_start; # this is inside Transcript obj
  }
  if ($label == $firstlabel) {
    return (1);
  }
  my ($start,$end,$strand)=($self->start(),$self->end(),$self->strand());
  my ($position,$in_pos,$out_pos,$coord_pos);
  my $length=$self->length;
  $coord_pos=$self->_inside_position($firstlabel);
  if ($self->valid($label)) { # if label is inside the Transcript
    $in_pos=$self->_inside_position($label);
    $position=$in_pos-$coord_pos+1;
    if ($position <= 0) {
      return ($position-1); # accounts for the missing of the 0 position
    }
  } else {
    if ($self->follows($end,$label)) { # label after end of transcript
      $out_pos=$self->{'seq'}->position($label,$end,$strand);
      #print "\n+++++++++DEBUG label $label FOLLOWS end $end outpos $out_pos coordpos $coord_pos";
      $position=$out_pos+$length-$coord_pos;
    } elsif ($self->follows($label,$start)) { # label before begin of transcript
      #print "\n+++++++++DEBUG label $label BEFORE start $start outpos $out_pos coordpos $coord_pos";
      $out_pos=$self->{'seq'}->position($label,$start,$strand);
      $position=$out_pos-$coord_pos+1;
    } else { # label is in intron (not valid, not after, not before)!
      $self->warn("Cannot give position of label pointing to intron according to CDS numbering!",1);
      return (0);
    }
  }
  return ($position);
}

sub seq {
  my $self=shift;
  my ($exon,$str);
  my @exons=$self->all_Exons();
  foreach $exon (@exons) {
    $str .= $exon->seq();
  }
  return $str;
}

sub length {
  my $self=shift;
  my ($exon,$length);
  my @exons=$self->all_Exons();
  foreach $exon (@exons) {
    $length += $exon->length();
  }
  return $length;
}

sub all_labels {
  my $self=shift;
  my ($exon,@labels);
  my @exons=$self->all_Exons();
  foreach $exon (@exons) {
    push (@labels,$exon->all_labels());
  }
  return @labels;
}

# redefined here so that it will retrieve effective subseq without introns
# otherwise it would have retrieved an underlying DNA (possibly with introns)
# subsequence
# Drawback: this is really bulky, label->position and then a call to
# subseq that will do the opposite position-> label
#
# one day this can be rewritten as the main one so that the normal subseq
# will rely on this one and hence avoid this double (useless and lengthy)
# conversion between labels and positions
sub old_labelsubseq {
  my ($self,$start,$length,$end)=@_;
  my ($pos1,$pos2);
  if ($start) {
    unless ($self->valid($start)) {
      $self->warn("Start label not valid"); return (-1);
    }
    $pos1=$self->position($start);
  }
  if ($end) {
    if ($end == $start) {
      $length=1;
    } else {
      unless ($self->valid($end)) {
	$self->warn("End label not valid"); return (-1);
      }
      unless ($self->follows($start,$end) == 1) {
	$self->warn("End label does not follow Start label!"); return (-1);
      }
      $pos2=$self->position($end);
      undef $length;
    }
  }
  return ($self->subseq($pos1,$pos2,$length));
}

# rewritten, eventually

sub labelsubseq {
  my ($self,$start,$length,$end,$unsecuremode)=@_;
  unless (defined $unsecuremode && 
	  $unsecuremode eq "unsecuremoderequested") 
  { # to skip security checks (faster)
    if ($start) {
      unless ($self->valid($start)) {
	$self->warn("Start label not valid"); return (-1);
      }
    } else {
      $start=$self->start;
    }
    if ($end) {
      if ($end == $start) {
	$length=1;
	undef $end;
      } else {
	undef $length; # end argument overrides length argument
	unless ($self->valid($end)) {
	  $self->warn("End label not valid"); return (-1);
	}
	unless ($self->follows($start,$end) == 1) {
	  $self->warn("End label does not follow Start label!"); return (-1);
	}
      }
    } else {
      $end=$self->end;
    }
  }
  my ($seq,$exon,$startexon,$endexon); my @exonlabels;
  my @exons=$self->all_Exons;
  EXONCHECK:
  foreach $exon (@exons) {
    if ((!(defined($startexon)))&&($exon->valid($start))) { # checks only if not yet found
      $startexon=$exon;
    }
    if ($exon->valid($end)) {
      $endexon=$exon;
    }
    if ((!(defined($seq)) && (defined($startexon)))) { # initializes only once
      if ((defined($endexon)) && ($endexon eq $startexon)) { # then perfect, we are finished
	if ($length) {
	  $seq = $startexon->labelsubseq($start,$length,undef,"unsecuremoderequested");


	  last EXONCHECK;
	} else {
	  $seq = $startexon->labelsubseq($start,undef,$end,"unsecuremoderequested");
	}
	last EXONCHECK;
      } else { # get up to the end of the exon
	$seq = $startexon->labelsubseq($start,undef,undef,"unsecuremoderequested");
      }
    }
    if (($startexon)&&($exon ne $startexon)) {
      if (defined($endexon)) { # we arrived to the last exon
	$seq .= $endexon->labelsubseq(undef,undef,$end,"unsecuremoderequested"); # get from the start of the exon
	last EXONCHECK;

      } elsif (defined($startexon)) { # we are in a whole-exon-in-the-middle case
	  $seq .= $exon->seq; # we add it completely to the seq
      } # else, we still have to reach the start point, exon useless, we move on
      if ($length) { # if length argument specified
	if (($seq && (CORE::length($seq) >= $length))) {
	  last EXONCHECK;
	}
      }
    }
  }
  if ($length) {
    return (substr($seq,0,$length));
  } else {
    return ($seq);
  }
}


# argument: label
# returns: the objref and progressive number of the Exon containing that label
# errorcode: -1
sub in_which_Exon {
  my ($self,$label)=@_;
  my ($count,$exon);
  my @exons=$self->all_Exons;
  foreach $exon (@exons) {
    $count++; # 1st exon is numbered "1"
    if ($exon->valid($label)) {
      return ($exon,$count)
    }
  }
  return (-1); # if nothing found
}

# recoded to exploit the new fast labelsubseq()
# valid only inside Transcript
sub subseq {
  my ($self,$pos1,$pos2,$length) = @_;
  my ($str,$startlabel,$endlabel);
  if (defined ($pos1)) {
    if ($pos1 == 0) {  # if position = 0 complain
      $self->warn("Position cannot be 0!"); return (-1);
    }
    if ((defined ($pos2))&&($pos1>$pos2)) {
      $self->warn("1st position($pos1) cannot be > 2nd position($pos2)!"); return (-1);
    }
    $startlabel=$self->label($pos1);
    unless ($self->valid($startlabel)) {
      $self->warn("Start label not valid"); return (-1);
    }
    if ($startlabel < 1) {
      $self->warn("position $pos1 not valid as start of subseq!"); return (-1);
    }
  } else {
    $startlabel=$self->start;
  }
  if (defined ($pos2)) {
    if ($pos2 == 0) {  # if position = 0 complain
      $self->warn("Position cannot be 0!"); return (-1);
    }
    undef $length;
    if ((defined ($pos1))&&($pos1>$pos2)) {
      $self->warn("1st position($pos1) cannot be > 2nd position($pos2)!"); return (-1);
    }
    $endlabel=$self->label($pos2);
    unless ($self->valid($endlabel)) {
      $self->warn("End label not valid"); return (-1);
    }
    if ($endlabel < 1) {
      $self->warn("position $pos2 not valid as end of subseq!"); return (-1);
    }
  } else {
    unless (defined ($length)) {
      $endlabel=$self->end;
    }
  }
  return ($self->labelsubseq($startlabel,$length,$endlabel,"unsecuremoderequested"));
}

# works only inside the transcript, complains if asked outside
sub old_subseq {
  my ($self,$pos1,$pos2,$length) = @_;
  my ($str,$startcount,$endcount,$seq,$seqlength);
  if (defined ($length)) {
    if ($length < 1) {
      $self->warn("No sense asking for a subseq of length < 1");
      return (-1);
    }
  }
  my $firstlabel=$self->coordinate_start; # this is inside Transcript obj
  my $coord_pos=$self->_inside_position($firstlabel); # TESTME old
  $seq=$self->seq;
  $seqlength=CORE::length($seq);
  unless (defined ($pos1)) {
    $startcount=1+$coord_pos-1; # i.e. coord_pos
  } else {
    if ($pos1 == 0) {  # if position = 0 complain
      $self->warn("Position cannot be 0!"); return (-1);
    } elsif ($pos1 < 0) {
      $pos1++;
    }
    if ((defined ($pos2))&&($pos1>$pos2)) {
      $self->warn("1st position ($pos1) cannot be > 2nd position ($pos2)!");
      return (-1);
    }
    $startcount=$pos1+$coord_pos-1;
  }
  unless (defined ($pos2)) {
  ;
  } else {
    if ($pos2 == 0) {  # if position = 0 complain
      $self->warn("Position cannot be 0!"); return (-1);
    } elsif ($pos2 < 0) {
      $pos2++;
    }
    if ((defined ($pos1))&&($pos1>$pos2)) {
      $self->warn("1st position ($pos1) cannot be > 2nd position ($pos2)!");
      return (-1);
    }
    $endcount=$pos2+$coord_pos-1;
    if ($endcount > $seqlength) {
      #print "\n###DEBUG###: pos1 $pos1 pos2 $pos2 coordpos $coord_pos endcount $endcount seqln $seqlength\n";
      $self->warn("Cannot access end position after the end of Transcript");
      return (-1);
    }
    $length=$endcount-$startcount+1;
  }
  #print "\n###DEBUG pos1 $pos1 pos2 $pos2 start $startcount end $endcount length $length coordpos $coord_pos\n";
  my $offset=$startcount-1;
  if ($offset < 0) {
    $self->warn("Cannot access startposition before the beginning of Transcript, returning from start",1); # ignorable
    return (substr($seq,0,$length));
  } elsif ($offset >= $seqlength) {
    $self->warn("Cannot access startposition after the end of Transcript");
    return (-1);
  } else {
    $str=substr($seq,$offset,$length);
    if (CORE::length($str) < $length) {
      $self->warn("Attention, cannot return the length requested ".
		  "for subseq",1) if $self->verbose > 0; # ignorable
    }
    return $str;
  }
}

# redefined so that it doesn't require other methods (after deletions) to
# reset it.
sub start {
  my $self = shift;
  my $exonsref=$self->{'exons'};
  my @exons=@{$exonsref};
  return ($exons[0]->start);
}

sub end {
  my $self = shift;
  my $exonsref=$self->{'exons'};
  my @exons=@{$exonsref};
  return ($exons[-1]->end);
}


# internal methods begin here

# returns: position of label in transcript's all_labels
#          with STARTlabel == 1
# errorcode 0 -> label not found
# argument: label
sub _inside_position {
  my ($self,$label)=@_;
  my ($start,$end,$strand)=($self->start(),$self->end(),$self->strand());
  my ($position,$checkme);
  my @labels=$self->all_labels;
  foreach $checkme (@labels) {
    $position++;
    if ($label == $checkme) {
      return ($position);
    }
  }
  return (0);
}

# returns 1 OK or 0 ERROR
# arguments: reference to array of Exon object references
sub _checkexons {
  my ($exon,$thisstart);
  my $self=$exon;
  my $exonsref=$_[0];
  my @exons=@{$exonsref};

  my $firstexon = $exons[0];

  unless (ref($firstexon) eq "Bio::LiveSeq::Exon") {
    $self->warn("Object not of class Exon");
    return (0);
  }
  my $strand = $firstexon->strand;

  my $prevend = $firstexon->end;
  shift @exons; # skip first one
  foreach $exon (@exons) {
    unless (ref($exon) eq "Bio::LiveSeq::Exon") { # object class check
      $self->warn("Object not of class Exon");
      return (0);
    }
    if ($exon->strand != $strand) { # strand consistency check
      $self->warn("Exons' strands not consistent when trying to create Transcript");
      return (0);
    }
    $thisstart = $exon->start;
    unless ($exon->{'seq'}->follows($prevend,$thisstart,$strand)) {
      $self->warn("Exons not in correct order when trying to create Transcript");
      return (0);
    }
    $prevend = $exon->end;
  }
  return (1);
}

=head2 get_Translation

  Title   : valid
  Usage   : $translation = $obj->get_Translation()
  Function: retrieves the reference to the object of class Translation (if any)
            attached to a LiveSeq object
  Returns : object reference
  Args    : none

=cut

sub get_Translation {
  my $self=shift;
  return ($self->{'translation'}); # this is set when Translation->new is called
}

# this checks so that deletion spanning multiple exons is
# handled accordingly and correctly
# arguments: begin and end label of a deletion
# this is called BEFORE any deletion in the chain
sub _deletecheck {
  my ($self,$startlabel,$endlabel)=@_;
  my $exonsref=$self->{'exons'};
  my @exons=@{$exonsref};
  my ($startexon,$endexon,$exon);
  $startexon=$endexon=0;
  foreach $exon (@exons) {
    if (($startexon == 0)&&($exon->valid($startlabel))) {
      $startexon=$exon; # exon containing start of deletion
    }
    if (($endexon == 0)&&($exon->valid($endlabel))) {
      $endexon=$exon; # exon containing end of deletion
    }
    if (($startexon)&&($endexon)) {
      last; # don't check further
    }
  }
  my $nextend=$self->label(2,$endlabel); # retrieve the next label
  my $prevstart=$self->label(-1,$startlabel); # retrieve the prev label

  if ($startexon eq $endexon) { # intra-exon deletion
    if (($startexon->start eq $startlabel) && ($startexon->end eq $endlabel)) {
      # let's delete the entire exon
      my @newexons;
      foreach $exon (@exons) {
	unless ($exon eq $startexon) {
	  push(@newexons,$exon);
	}
      }
      $self->{'exons'}=\@newexons;
    } elsif ($startexon->start eq $startlabel) { # special cases
      $startexon->{'start'}=$nextend; # set a new start of exon
    } elsif ($startexon->end eq $endlabel) {
      $startexon->{'end'}=$prevstart; # set a new end of exon
    } else {
      return; # no problem
    }
  } else { # two new exons to be created, inter-exons deletion
    my @newexons;
    my $exonobj;
    my $dna=$self->{'seq'};
    my $strand=$self->strand;
    my $notmiddle=1; # flag for skipping exons in the middle of deletion
    foreach $exon (@exons) {
      if ($exon eq $startexon) {
	$exonobj=Bio::LiveSeq::Exon->new('-seq'=>$dna,'-start'=>$exon->start,'-end'=>$prevstart,'-strand'=>$strand); # new partial exon
	push(@newexons,$exonobj);
	$notmiddle=0; # now we enter totally deleted exons
      } elsif ($exon eq $endexon) {
	$exonobj=Bio::LiveSeq::Exon->new('-seq'=>$dna,'-start'=>$nextend,'-end'=>$exon->end,'-strand'=>$strand); # new partial exon
	push(@newexons,$exonobj);
	$notmiddle=1; # exiting totally deleted exons
      } else {
	if ($notmiddle) { # if before or after exons with deletion
	  push(@newexons,$exon); 
	}# else skip them
      }
    }
    $self->{'exons'}=\@newexons;
  }
}

=head2 translation_table

 Title   : translation_table
 Usage   : $name = $obj->translation_table;
         : $name = $obj->translation_table(11);
 Function: Returns or sets the translation_table used for translating the
           transcript.
           If it has never been set, it will return undef.
 Returns : an integer

=cut

sub translation_table {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'translation_table'} = $value;
  }
  unless (exists $self->{'translation_table'}) {
    return;
  } else {
    return $self->{'translation_table'};
  }
}

=head2 frame

 Title   : frame
 Usage   : $frame = $transcript->frame($label);
 Function: Returns the frame of a particular nucleotide.
           Frame can be 0 1 or 2 and means the position in the codon triplet
           of the particulat nucleotide. 0 is the first codon_position.
           Codon_position (1 2 3) is simply frame+1.
           If the label asked for is not inside the Transcript, -1 will be
           returned.
 Args    : a label
 Returns : 0 1 or 2
 Errorcode -1

=cut

# args: label
# returns: frame of nucleotide (0 1 2)
# errorcode: -1
sub frame {
  my ($self,$inputlabel)=@_;
  my @labels=$self->all_labels;
  my ($label,$frame,$count);
  foreach $label (@labels) {
    if ($inputlabel == $label) {
      return ($count % 3);
    }
    $count++; # 0 1 2 3 4....
  }
  return (-1); # label not found amid Transcript labels
}

1;
