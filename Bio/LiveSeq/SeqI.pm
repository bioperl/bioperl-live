#
# bioperl module for Bio::LiveSeq::SeqI
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

Bio::LiveSeq::SeqI - Abstract sequence interface class for LiveSeq

=head1 SYNOPSIS

  # documentation needed

=head1 DESCRIPTION

This class implements BioPerl PrimarySeqI interface for Live Seq objects.

One of the main difference in LiveSequence compared to traditional
"string" sequences is that coordinate systems are flexible. Typically
gene nucleotide numbering starts from 1 at the first character of the
initiator codon (A in ATG). This means that negative positions are
possible and common!

Secondly, the sequence manipulation methods do not return a new
sequence object but change the current object. The current status can
be written out to BioPerl sequence objects.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

Some note on the terminology/notation of method names:
 label: a unique pointer to a single nucleotide
 position: the position of a nucleotide according to a particular coordinate
           system (e.g. counting downstream from a particular label taken as
           number 1)
 base: the one letter code for a nucleotide (i.e.: "a" "t" "c" "g")

       a base is the "value" that an "element" of a "chain" can assume
         (see documentation on the Chain datastructure if interested)

=cut

#'
# Let the code begin...

package Bio::LiveSeq::SeqI;
use strict;
use Bio::Tools::CodonTable; # for the translate() function

use base qw(Bio::Root::Root Bio::LiveSeq::ChainI Bio::PrimarySeqI);

=head2 seq

 Title   : seq
 Usage   : $string    = $obj->seq()
 Function: Returns the complete sequence of an object as a string of letters.
           Suggested cases are upper case for proteins and lower case for
           DNA sequence (IUPAC standard),
 Returns : a string


=cut

sub seq {
  my $self = shift;
  my ($start,$end) = ($self->start(),$self->end());
  if ($self->strand() == 1) {
    return $self->{'seq'}->down_chain2string($start,undef,$end);
  } else { # reverse strand
    my $str = $self->{'seq'}->up_chain2string($start,undef,$end);
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    return $str;
  }
}

=head2 all_labels

 Title   : all_labels
 Usage   : @labels = $obj->all_labels()
 Function: all the labels of every nucleotide an object is composed of
 Returns : an array of labels
 Args    : none

=cut

sub all_labels {
  my $self = shift;
  my ($start,$end) = ($self->start(),$self->end());
  my $labels;
  if ($self->strand() == 1) {
    $labels=$self->{'seq'}->down_labels($start,$end);
  } else {
    $labels=$self->{'seq'}->up_labels($start,$end);
  }
  return (@{$labels});
}

=head2 labelsubseq

  Title   : labelsubseq
  Usage   : $dna->labelsubseq();
          : $dna->labelsubseq($startlabel);
          : $dna->labelsubseq($startlabel,$length);
          : $dna->labelsubseq($startlabel,undef,$endlabel);
  e.g.    : $dna->labelsubseq(4,undef,8);
  Function: prints the sequence as string. The difference between labelsubseq
            and normal subseq is that it uses /labels/ as arguments, instead
            than positions. This allows for faster and more efficient lookup,
            skipping the (usually) lengthy conversion of positions into labels.
            This is especially useful for manipulating with high power
            LiveSeq objects, knowing the labels and exploiting their
            usefulness.
  Returns : a string
  Errorcode -1
  Args    : without arguments it returns the entire sequence
            with a startlabel it returns the sequence downstream that label
            if a length is specified, it returns only that number of bases
            if an endlabel is specified, it overrides the length argument
             and prints instead up to that label (included)
  Defaults: $startlabel defaults to the beginning of the entire sequence
            $endlabel defaults to the end of the entire sequence

=cut

# NOTE: unsecuremode is to be used /ONLY/ if sure of the start and end labels, especially that they follow each other in the correct order!!!!

sub labelsubseq {
  my ($self,$start,$length,$end,$unsecuremode) = @_;
  if (defined $unsecuremode && $unsecuremode eq "unsecuremoderequested")
  { # to skip security checks (faster)
    unless ($start) {
      $start=$self->start;
    }
    if ($end) {
      if ($end == $start) {
	$length=1;
	undef $end;
      } else {
	undef $length;
      }
    } else {
      unless ($length) {
	$end=$self->end;
      }
    }
  } else {
    if ($start) {
      unless ($self->{'seq'}->valid($start)) {
	$self->warn("Start label not valid"); return (-1);
      }
    }
    if ($end) {
      if ($end == $start) {
	$length=1;
	undef $end;
      } else {
	unless ($self->{'seq'}->valid($end)) {
	  $self->warn("End label not valid"); return (-1);
	}
	unless ($self->follows($start,$end) == 1) {
	  $self->warn("End label does not follow Start label!"); return (-1);
	}
	undef $length;
      }
    }
  }
  if ($self->strand() == 1) {
    return $self->{'seq'}->down_chain2string($start,$length,$end);
  } else { # reverse strand
    my $str = $self->{'seq'}->up_chain2string($start,$length,$end);
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    return $str;
  }
}

=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10,40);
         : $substring = $obj->subseq(10,undef,4);
 Function: returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence

           Start cannot be larger than end but can be equal.

           Allows for negative numbers $obj->subseq(-10,-1). By
           definition, there is no 0!
                       -5  -1 1   5
                gctagcgcccaac atggctcgctg

           This allows one to retrieve sequences upstream from given position.

           The precedence is from left to right: if END is given LENGTH is
           ignored.

 Examples: $obj->subseq(-10,undef,10) returns 10 elements before position 1
           $obj->subseq(4,8) returns elements from the 4th to the 8th, inclusive

 Returns : a string
 Errorcode: -1
 Args    : start,  integer, defaults to start of the sequence
           end,    integer, '' or undef, defaults to end of the sequence
           length, integer, '' or undef
           an optional strand (1 or -1) 4th argument
            if strand argument is not given, it will default to the object
            argment. This argument is useful when a call is issued from a child
            of a parent object containing the subseq method

=cut

#'
# check the fact about reverse strand!
# is it feasible? Is it correct? Should we do it? How about exons? Does it
# work when you ask subseq of an exon?
# eliminated now (Mon night)
sub subseq {
  ##my ($self,$pos1,$pos2,$length,$strand) = @_;
  my ($self,$pos1,$pos2,$length,$strand) = @_;
  ##unless (defined ($strand)) { # if optional [strand] argument not given
  ##  $strand=$self->strand;
  ##}
  $strand=$self->strand;
  my ($str,$startlabel,$endlabel);
  if (defined ($length)) {
    if ($length < 1) {
      $self->warn("No sense asking for a subseq of length < 1");
      return (-1);
    }
  }
  unless (defined ($pos1)) {
    #print "\n##### DEBUG pos1 not defined\n";
    $startlabel=$self->start;
  } else {
    if ($pos1 == 0) {  # if position = 0 complain
      $self->warn("Position cannot be 0!"); return (-1);
    }
    ##if ($strand == 1) { # CHECK THIS!
      if ((defined ($pos2))&&($pos1>$pos2)) {
	$self->warn("1st position($pos1) cannot be > 2nd position($pos2)!"); return (-1);
      }
    ##} else { # CHECK THIS!
    ##  if ((defined ($pos2))&&($pos1<$pos2)) {
##	$self->warn("1st position($pos1) cannot be < 2nd position($pos2) on reverse strand!)"; return (-1);
    ##  }
    ##}
    $startlabel=$self->label($pos1);
    if ($startlabel < 1) {
      $self->warn("position $pos1 not valid as start of subseq!"); return (-1);
    }
  }
  unless (defined ($pos2)) {
    #print "\n##### pos2 not defined\n";
    unless (defined ($length)) {
      $endlabel=$self->end;
    }
  } else {
    if ($pos2 == 0) {  # if position = 0 complain
      $self->warn("Position cannot be 0!"); return (-1);
    }
    undef $length;
    ##if ($strand == 1) { # CHECK THIS!
      if ((defined ($pos1))&&($pos1>$pos2)) {
	$self->warn("1st position($pos1) cannot be > 2nd position($pos2)!"); return (-1);
      }
    ##} else { # CHECK THIS!
    ##  if ((defined ($pos1))&&($pos1<$pos2)) {
##	$self->warn("1st position($pos1) cannot be < 2nd position($pos2) on reverse strand!"); return (-1);
    ##  }
    ##}
    $endlabel=$self->label($pos2);
    if ($endlabel < 1) {
      $self->warn("position $pos2 not valid as end of subseq!"); return (-1);
    }
  }
  #print "\n    ####DEBUG: start $startlabel end $endlabel length $length strand $strand\n";

  if ($strand == 1) {
    $str = $self->{'seq'}->down_chain2string($startlabel,$length,$endlabel);
  } else { # reverse strand
    $str = $self->{'seq'}->up_chain2string($startlabel,$length,$endlabel);
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  }
  return $str;
}

=head2 length

  Title   : length
  Usage   : $seq->length();
  Function: returns the number of nucleotides (or the number of aminoacids)
            in the entire sequence
  Returns : an integer
  Errorcode -1
  Args    : none

=cut

sub length {
  my $self=shift;
  my ($start,$end,$strand)=($self->start(),$self->end(),$self->strand());
  if ($strand == 1) {
    return $self->{'seq'}->down_subchain_length($start,$end);
  } else {
    return $self->{'seq'}->up_subchain_length($start,$end);
  }
}

=head2 display_id

 Title   : display_id
 Usage   : $id_string = $obj->display_id();
 Function: returns the display id, alias the common name of the object

           The semantics of this is that it is the most likely string
           to be used as an identifier of the sequence, and likely to
           have "human" readability.  The id is equivalent to the ID
           field of the GenBank/EMBL databanks and the id field of the
           Swissprot/sptrembl database. In fasta format, the >(\S+) is
           presumed to be the id, though some people overload the id
           to embed other information.

 See also: accession_number
 Returns : a string
 Args    : none

=cut

sub display_id {
  my ($self,$value) = @_;
  if(defined $value) {
    $self->{'display_id'} = $value;
  }
  return $self->{'display_id'};
}


=head2 accession_number

 Title   : accession_number
 Usage   : $unique_biological_key = $obj->accession_number;
 Function: Returns the unique biological id for a sequence, commonly
           called the accession_number.
           Notice that primary_id() provides the unique id for the
           implemetation, allowing multiple objects to have the same accession
           number in a particular implementation.

           For objects with no accession_number this method returns "unknown".
 Returns : a string
 Args    : none

=cut

sub accession_number {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'accession_number'} = $value;
  }
  unless (exists $self->{'accession_number'}) {
    return "unknown";
  } else {
    return $self->{'accession_number'};
  }
}

=head2 primary_id

 Title   : primary_id
 Usage   : $unique_implementation_key = $obj->primary_id;
 Function: Returns the unique id for this object in this
           implementation. This allows implementations to manage their own
           object ids in a way the implementation can control. Clients can
           expect one id to map to one object.

           For sequences with no primary_id, this method returns
           a stringified memory location.

 Returns : A string
 Args    : None

=cut


sub primary_id {
  my ($self,$value) = @_;
  if(defined $value) {
    $self->{'primary_id'} = $value;
  }
  unless (exists $self->{'primary_id'}) {
    return "$self";
  } else {
    return $self->{'primary_id'};
  }
}

=head2 change

 Title   : change
 Usage   : $substring = $obj->change('AA', 10);
 Function: changes, modifies, mutates the LiveSequence
 Examples:
        $obj->change('',   10);      delete nucleotide #10
        $obj->change('',   10, 2);   delete two nucleotides starting from #10
        $obj->change('G',  10);      change nuc #10 to 'G'
        $obj->change('GA', 10, 4);   replace #10 and 3 following with 'GA'
        $obj->change('GA', 10, 2));  is same as $obj->change('GA',  10);
        $obj->change('GA', 10, 0 );  insert 'GA' before nucleotide at #10
        $obj->change('GA', 10, 1);   GA inserted before #10, #10 deleted
        $obj->change('GATC', 10, 2); GATC inserted before #10, #10&#11 deleted
        $obj->change('GATC', 10, 6); GATC inserted before #10, #10-#15 deleted


 Returns : a string of deleted bases (if any) or 1 (everything OK)
 Errorcode: -1
 Args    : seq,    string, or '' ('' = undef = 0 = deletion)
           start,  integer
           length, integer (optional)

=cut

sub change {
  &positionchange;
}

=head2 positionchange

 Title   : positionchange
 Function: Exactly like change. I.e. change() defaults to positionchange()

=cut

sub positionchange {
  my ($self,$newseq,$position,$length)=@_;
  unless ($position) {
    $self->warn("Position not given or position 0");
    return (-1);
  }
  my $label=$self->label($position);
  unless ($label > 0) { # label not found or error
    $self->warn("No valid label found at that position!");
    return (-1);
  }
  return ($self->labelchange($newseq,$label,$length));
}

=head2 labelchange

 Title   : labelchange
 Function: Exactly like change but uses a /label/ instead than a position
           as second argument. This allows for multiple changes in a LiveSeq
           without the burden of recomputing positions. I.e. for a multiple
           change in two different points of the LiveSeq, the approach would
           be the following: fetch the correct labels out of the two different
           positions (method: label($position)) and then use the labelchange()
           method to modify the sequence using those labels instead than
           relying on the positions (that would have modified after the
           first change).

=cut

sub labelchange {
  my ($self,$newseq,$label,$length)=@_;
  unless ($self->valid($label)) {
    if ($self->{'seq'}->valid($label)) {
       #$self->warn("Label \'$label\' not valid for executing a LiveSeq change for the object asked but it's ok for DNAlevel change, reverting to that");
      shift @_;
      return($self->{'seq'}->labelchange(@_));
    } else {
      $self->warn("Label \'$label\' not valid for executing a LiveSeq change");
      return (-1);
    }
  }
  unless ($newseq) { # it means this is a simple deletion
    if (defined($length)) {
      unless ($length >= 0) {
	$self->warn("No sense having length < 0 in a deletion");
	return (-1);
      }
    } else {
      $self->warn("Length not defined for deletion!");
      return (-1);
    }
    return $self->_delete($label,$length);
  }
  my $newseqlength=CORE::length($newseq);
  if (defined($length)) {
    unless ($length >= 0) {
      $self->warn("No sense having length < 0 in a change()");
      return (-1);
    }
  } else {
    $length=$newseqlength; # defaults to pointmutation(s)
  }
  if ($length == 0) { # it means this is a simple insertion, length def&==0
    my ($insertbegin,$insertend)=$self->_praeinsert($label,$newseq);
    if ($insertbegin == -1) {
      return (-1);
    } else {
      return (1);
    }
  }
  if ($newseqlength == $length) { # it means this is simple pointmutation(s)
    return $self->_mutate($label,$newseq,$length);
  }
  # if we arrived here then change is complex mixture
  my $strand=$self->strand();
  my $afterendlabel=$self->label($length+1,$label,$strand); # get the label at $length+1 positions after $label
  unless ($afterendlabel > 0) { # label not found or error
    $self->warn("No valid afterendlabel found for executing the complex mutation!");
    return (-1);
  }
  my $deleted=$self->_delete($label,$length); # first delete length nucs
  if ($deleted eq -1) { # if errors
    return (-1);
  } else { # then insert the newsequence
    my ($insertbegin,$insertend)=$self->_praeinsert($afterendlabel,$newseq);
    if ($insertbegin == -1) {
      return (-1);
    } else {
      return (1);
    }
  }
}

# internal methods for change()

# arguments: label for beginning of deletion, new sequence to insert
# returns: labels of beginning and end of the inserted sequence
# errorcode: -1
sub _praeinsert {
  my ($self,$label,$newseq)=@_;
  my ($insertbegin,$insertend);
  my $strand=$self->strand();
  if ($strand == 1) {
    ($insertbegin,$insertend)=($self->{'seq'}->praeinsert_string($newseq,$label));
  } else { # since it's reverse strand and we insert in forward direction....
    $newseq=reverse($newseq);
    $newseq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/; # since it's reverse strand we get the complementary bases
    ($insertend,$insertbegin)=($self->{'seq'}->postinsert_string($newseq,$label));
  }
  if (($insertbegin==0)||($insertend==0)) {
    $self->warn("Some error occurred while inserting!");
    return (-1);
  } else {
    return ($insertbegin,$insertend);
  }
}

# arguments: label for beginning of deletion, length of deletion
# returns: string of deleted bases
# errorcode: -1
sub _delete {
  my ($self,$label,$length)=@_;
  my $strand=$self->strand();
  my $endlabel=$self->label($length,$label,$strand); # get the label at $length positions after $label
  unless ($endlabel > 0) { # label not found or error
    $self->warn("No valid endlabel found for executing the deletion!");
    return (-1);
  }
  # this is important in Transcript to fix exon structure
  $self->_deletecheck($label,$endlabel);
  my $deletedseq;
  if ($strand == 1) {
    $deletedseq=$self->{'seq'}->splice_chain($label,undef,$endlabel);
  } else {
    $deletedseq=$self->{'seq'}->splice_chain($endlabel,undef,$label);
    $deletedseq=reverse($deletedseq); # because we are on reverse strand and we cut anyway
                         # in forward direction
    $deletedseq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/; # since it's reverse strand we get the complementary bases
  }
  return ($deletedseq);
}

# empty function, overridden in Transcript, not useful here
sub _deletecheck {
}

# arguments: label for beginning of mutation, newsequence, number of mutations
# returns: 1 all OK
# errorcode: -1
sub _mutate {
  my ($self,$label,$newseq,$length)=@_; # length is equal to length(newseq)
  my ($i,$base,$nextlabel);
  my @labels; # array of labels
  my $strand=$self->strand();
  if ($length == 1) { # special cases first
    @labels=($label);
  } else {
    my $endlabel=$self->label($length,$label,$strand); # get the label at $length positions after $label
    unless ($endlabel > 0) { # label not found or error
      $self->warn("No valid endlabel found for executing the mutation!");
      return (-1);
    }
    if ($length == 2) { # another special case
      @labels=($label,$endlabel);
    } else { # more than 3 bases changed
      # this wouldn't work for Transcript
      #my $labelsarrayref;
      #if ($strand == 1) {
	#$labelsarrayref=$self->{'seq'}->down_labels($label,$endlabel);
      #} else {
	#$labelsarrayref=$self->{'seq'}->up_labels($label,$endlabel);
      #}
      #@labels=@{$labelsarrayref};
      #if ($length != scalar(@labels)) { # not enough labels returned
	#$self->warn("Not enough valid labels found for executing the mutation!");
	#return (-1);
      #}

      # this should be more general
      @labels=($label); # put the first one
      while ($label != $endlabel) {
	$nextlabel=$self->label(2,$label,$strand); # retrieve the next label
	push (@labels,$nextlabel);
	$label=$nextlabel; # move on reference
      }
    }
  }
  if ($strand == -1) { # only for reverse strand
    $newseq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/; # since it's reverse strand we get the complementary bases
  }
  my $errorcheck; # if not equal to $length after summing for all changes, error did occurr
  $i = 0;
  foreach $base (split(//,$newseq)) {
    $errorcheck += $self->{'seq'}->set_value_at_label($base,$labels[$i]);
    $i++;
  }
  if ($errorcheck != $length) {
    $self->warn("Some error occurred while mutating!");
    return (-1);
  } else {
    return (1);
  }
}

=head2 valid

  Title   : valid
  Usage   : $boolean = $obj->valid($label)
  Function: tests if a label exists inside the object
  Returns : boolean
  Args    : label

=cut

# argument: label
# returns: 1 YES 0 NO
sub valid {
  my ($self,$label)=@_;
  my $checkme;
  my @labels=$self->all_labels;
  foreach $checkme (@labels) {
    if ($label == $checkme) {
      return (1); # found
    }
  }
  return (0); # not found
}


=head2 start

  Title   : start
  Usage   : $startlabel=$obj->start()
  Function: returns the label of the first nucleotide of the object (exon, CDS)
  Returns : label
  Args    : none

=cut

sub start {
  my ($self) = @_;
  return $self->{'start'}; # common for all classes BUT DNA (which redefines it) and Transcript (that takes the information from the Exons)
}

=head2 end

  Title   : end
  Usage   : $endlabel=$obj->end()
  Function: returns the label of the last nucleotide of the object (exon, CDS)
  Returns : label
  Args    : none

=cut

sub end {
  my ($self) = @_;
  return $self->{'end'};
}

=head2 strand

  Title   : strand
  Usage   : $strand=$obj->strand()
            $obj->strand($strand)
  Function: gets or sets strand information, being 1 or -1 (forward or reverse)
  Returns : -1 or 1
  Args    : none OR -1 or 1

=cut

sub strand {
  my ($self,$strand) = @_;
  if ($strand) {
    if (($strand != 1)&&($strand != -1)) {
      $self->warn("strand information not changed because strand identifier not valid");
    } else {
      $self->{'strand'} = $strand;
    }
  }
  return $self->{'strand'};
}

=head2 alphabet

 Title   : alphabet
 Usage   : if( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence being one of
           'dna', 'rna' or 'protein'. This is case sensitive.

 Returns : a string either 'dna','rna','protein'.
 Args    : none

=cut


sub alphabet {
  my %valid_type = map {$_, 1} qw( dna rna protein );
  my ($self,$value) = @_;
  if (defined $value) {
    $value = 'dna' if $value =~ /dna/i;
    $value = 'rna' if $value =~ /rna/i;
    unless ( $valid_type{$value} ) {
      $self->warn("Molecular type '$value' is not a valid type");
    }
    $self->{'alphabet'} = $value;
  }
  return $self->{'alphabet'};
}

=head2 coordinate_start

  Title   : coordinate_start
  Usage   : $coordstartlabel=$obj->coordinate_start()
          : $coordstartlabel=$obj->coordinate_start($label)
  Function: returns and optionally sets the first label of the coordinate
            system used
            For some objects only labels inside the object or in frame (for
            Translation objects) will be allowed to get set as coordinate start

  Returns : label. It returns 0 if label not found.
  Errorcode -1
  Args    : an optional reference $label that is position 1

=cut


sub coordinate_start {
  my ($self,$label) = @_;
  if ($label) {
    if ($self->valid($label)) {
      $self->{'coordinate_start'} = $label;
    } else {
      $self->warn("The label you are trying to set as coordinate_start is not valid for this object");
    }
  }
  my $coord_start = $self->{'coordinate_start'};
  if ($coord_start) {
    return $coord_start;
  } else {
    return $self->start();
  }
}

=head2 label

  Title   : label
  Usage   : $seq->label($position)
          : $seq->label($position,$firstlabel)
  Examples: $nextlabel=$seq->label(2,$label) -> retrieves the following label
          : $prevlabel=$seq->label(-1,$label) -> retrieves the preceding label

  Function: returns the label of the nucleotide at $position from current
            coordinate start
  Returns : a label. It returns 0 if label not found.
  Errorcode -1
  Args    : a position,
            an optional reference $firstlabel that is to be used as position 1
            an optional strand (1 or -1) argument
             if strand argument is not given, it will default to the object
             argument. This argument is useful when a call is issued from a child
             of a parent object containing the subseq method

=cut


sub label {
  my ($self,$position,$firstlabel,$strand)=@_;
  my $label;
  unless (defined ($firstlabel)) {
    $firstlabel=$self->coordinate_start;
  }
  unless ($position) {  # if position = 0 complain ?
    $self->warn("Position not given or position 0");
    return (-1);
  }
  unless (defined ($strand)) { # if optional [strand] argument not given
    $strand=$self->strand;
  }
  if ($strand == 1) {
    if ($position > 0) {
      $label=$self->{'seq'}->down_get_label_at_pos($position,$firstlabel)
    } else { # if < 0
      $label=$self->{'seq'}->up_get_label_at_pos(1 - $position,$firstlabel)
    }
  } else {
    if ($position > 0) {
      $label=$self->{'seq'}->up_get_label_at_pos($position,$firstlabel)
    } else { # if < 0
      $label=$self->{'seq'}->down_get_label_at_pos(1 - $position,$firstlabel)
    }
  }
  return $label;
}


=head2 position

  Title   : position
  Usage   : $seq->position($label)
          : $seq->position($label,$firstlabel)
  Function: returns the position of nucleotide at $label
  Returns : the position of the label from current coordinate start
  Errorcode 0
  Args    : a label pointing to a certain nucleotide (e.g. start of exon)
            an optional "firstlabel" as reference to count from
            an optional strand (1 or -1) argument
             if strand argument is not given, it will default to the object
             argument. This argument is useful when a call is issued from a child
             of a parent object containing the subseq method

=cut


sub position {
  my ($self,$label,$firstlabel,$strand)=@_;
  unless (defined ($strand)) { # if optional [strand] argument not given
    $strand=$self->strand;
  }
  unless (defined ($firstlabel)) {
    $firstlabel=$self->coordinate_start;
  }
  unless ($self->valid($label)) {
    $self->warn("label not valid");
    return (0);
  }
  if ($firstlabel == $label) {
    return (1);
  }
  my ($coordpos,$position0,$position);
  $position0=$self->{'seq'}->down_get_pos_of_label($label);
  $coordpos=$self->{'seq'}->down_get_pos_of_label($firstlabel);
  $position=$position0-$coordpos+1;
  if ($position <= 0) {
    $position--;
  }
  if ($strand == -1) {
    #print "\n----------DEBUGSEQPOS label $label firstlabel $firstlabel strand $strand: position=",1-$position;
    return (1-$position);
  } else {
    #print "\n----------DEBUGSEQPOS label $label firstlabel $firstlabel strand $strand: position=",$position;
    return ($position);
  }
}

=head2 follows

  Title   : follows
  Usage   : $seq->follows($firstlabel,$secondlabel)
          : $seq->follows($firstlabel,$secondlabel,$strand)
  Function: checks if SECONDlabel follows FIRSTlabel, undependent of the strand
            i.e. it checks downstream for forward strand and
            upstream for reverse strand
  Returns : 1 or 0
  Errorcode -1
  Args    : two labels
            an optional strand (1 or -1) argument
             if strand argument is not given, it will default to the object
             argument. This argument is useful when a call is issued from a child
             of a parent object containing the subseq method

=cut

#'
# wraparound to is_downstream and is_upstream that chooses the correct one
# depending on the strand
sub follows {
  my ($self,$firstlabel,$secondlabel,$strand)=@_;
  unless (defined ($strand)) { # if optional [strand] argument not given
    $strand=$self->strand;
  }
  if ($strand == 1) {
    return ($self->{'seq'}->is_downstream($firstlabel,$secondlabel));
  } else {
    return ($self->{'seq'}->is_upstream($firstlabel,$secondlabel));
  }
}
#
#=head2 translate
#
# Title   : translate
# Usage   : $protein_seq = $obj->translate
# Function: Provides the translation of the DNA sequence
#	    using full IUPAC ambiguities in DNA/RNA and amino acid codes.
#
#	    The resulting translation is identical to EMBL/TREMBL database
#	    translations.
#
# Returns : a string
# Args    : character for terminator (optional) defaults to '*'
#	    character for unknown amino acid (optional) defaults to 'X'
#	    frame (optional) valid values 0, 1, 3, defaults to 0
#	    codon table id (optional) defaults to 1
#
#=cut
#
#sub translate {
#  my ($self) = shift;
#  return ($self->translate_string($self->seq,@_));
#}
#
#=head2 translate_string
#
# Title   : translate_string
# Usage   : $protein_seq = $obj->translate_string("attcgtgttgatcgatta");
# Function: Like translate, but can be used to translate subsequences after
#	    having retrieved them as string.
# Args    : 1st argument is a string. Optional following arguments: like in
#	    the translate method
#
#=cut
#
#
#sub translate_string {
#  my($self) = shift;
#  my($seq) = shift;
#  my($stop, $unknown, $frame, $tableid) = @_;
#  my($i, $len, $output) = (0,0,'');
#  my($codon)   = "";
#  my $aa;
#
#
#  ## User can pass in symbol for stop and unknown codons
#  unless(defined($stop) and $stop ne '')    { $stop = "*"; }
#  unless(defined($unknown) and $unknown ne '') { $unknown = "X"; }
#  unless(defined($frame) and $frame ne '') { $frame = 0; }
#
#  ## the codon table ID
#  if ($self->translation_table) {
#    $tableid = $self->translation_table;
#  }
#  unless(defined($tableid) and $tableid ne '')    { $tableid = 1; }
#
#  ##Error if monomer is "Amino"
#  $self->warn("Can't translate an amino acid sequence.")
#      if (defined $self->alphabet && $self->alphabet eq 'protein');
#
#  ##Error if frame is not 0, 1 or 2
#  $self->warn("Valid values for frame are 0, 1, 2, not [$frame].")
#      unless ($frame == 0 or $frame == 1 or $frame == 2);
#
#  #thows a warning if ID is invalid
#  my $codonTable = Bio::Tools::CodonTable->new( -id => $tableid);
#
#  # deal with frame offset.
#  if( $frame ) {
#      $seq = substr ($seq,$frame);
#  }
#
#  for $codon ( grep { CORE::length == 3 } split(/(.{3})/, $seq) ) {
#      my $aa = $codonTable->translate($codon);
#      if ($aa eq '*') {
#	    $output .= $stop;
#      }
#      elsif ($aa eq 'X') {
#	    $output .= $unknown;
#      }
#      else {
#	   $output .= $aa ;
#      }
#  }
#  #if( substr($output,-1,1) eq $stop ) {
#  #    chop $output;
#  #}
#
#  return ($output);
#}

=head2 gene

 Title   : gene
 Usage   : my $gene=$obj->gene;
 Function: Gets or sets the reference to the LiveSeq::Gene object.
           Objects that are features of a LiveSeq Gene will have this
           attribute set automatically.

 Returns : reference to an object of class Gene
 Note    : if Gene object is not set, this method will return 0;
 Args    : none or reference to object of class Bio::LiveSeq::Gene

=cut

sub gene {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'gene'} = $value;
  }
  unless (exists $self->{'gene'}) {
    return (0);
  } else {
    return $self->{'gene'};
  }
}

=head2 obj_valid

 Title   : obj_valid
 Usage   : if ($obj->obj_valid) {do something;}
 Function: Checks if start and end labels are still valid for the ojbect,
           i.e. tests if the LiveSeq object is still valid
 Returns : boolean
 Args    : none

=cut

sub obj_valid {
  my $self=shift;
  unless (($self->{'seq'}->valid($self->start()))&&($self->{'seq'}->valid($self->end()))) {
    return (0);
  }
  return (1);
}

=head2 name

 Title   : name
 Usage   : $name = $obj->name;
         : $name = $obj->name("ABCD");
 Function: Returns or sets the name of the object.
           If there is no name, it will return "unknown";
 Returns : A string
 Args    : None

=cut

sub name {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'name'} = $value;
  }
  unless (exists $self->{'name'}) {
    return "unknown";
  } else {
    return $self->{'name'};
  }
}

=head2 desc

 Title   : desc
 Usage   : $desc = $obj->desc;
         : $desc = $obj->desc("ABCD");
 Function: Returns or sets the description of the object.
           If there is no description, it will return "unknown";
 Returns : A string
 Args    : None

=cut

sub desc {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'desc'} = $value;
  }
  unless (exists $self->{'desc'}) {
    return "unknown";
  } else {
    return $self->{'desc'};
  }
}

=head2 source

 Title   : source
 Usage   : $name = $obj->source;
         : $name = $obj->source("Homo sapiens");
 Function: Returns or sets the organism that is source of the object.
           If there is no source, it will return "unknown";
 Returns : A string
 Args    : None

=cut

sub source {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'source'} = $value;
  }
  unless (exists $self->{'source'}) {
    return "unknown";
  } else {
    return $self->{'source'};
  }
}

sub delete_Obj {
  my $self = shift;
  my @values= values %{$self};
  my @keys= keys %{$self};

  foreach my $key ( @keys ) {
    delete $self->{$key};
  }
  foreach my $value ( @values ) {
    if (index(ref($value),"LiveSeq") != -1) { # object case
      eval {
	# delete $self->{$value};
	$value->delete_Obj;
      };
    } elsif (index(ref($value),"ARRAY") != -1) { # array case
      my @array=@{$value};
      my $element;
      foreach $element (@array) {
	eval {
	  $element->delete_Obj;
	};
      }
    } elsif (index(ref($value),"HASH") != -1) { # object case
      my %hash=%{$value};
      my $element;
      foreach $element (%hash) {
	eval {
	  $element->delete_Obj;
	};
      }
    }
  }
  return(1);
}

1;
