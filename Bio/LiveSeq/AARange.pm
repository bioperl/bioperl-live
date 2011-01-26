#
# bioperl module for Bio::LiveSeq::AARange
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

Bio::LiveSeq::AARange - AARange abstract class for LiveSeq

=head1 SYNOPSIS

  #documentation needed

=head1 DESCRIPTION

This is used as possible parent for aminoacid range object classes.
Or it can be used straight away to define aminoacid ranges.  The idea
is that the ranges defined are attached to a Translation object and
they refer to its coordinate-system when they are first created (via
the new() method).  When they are created they are anyway linked to
the underlying DNA LiveSeq by way of the LiveSeq labels. This allows
to preserve the ranges even if the numbering changes in the
Translation due to deletions or insertions.

The protein sequence associated with the AARange can be accessed via
the usual seq() or subseq() methods.

The start and end of the AARange in protein coordinate system can be
fetched with aa_start() and aa_end() methods. Note: the behaviour of
these methods would be influenced by the coordinate_start set in the
corresponding Translation object. This can be desirable but can also
lead to confusion if the coordinate_start had been changed and the
original position of the AARange was to be retrieved.

start() and end() methods of the AARange will point to the labels
identifying the first nucleotide of the first and last triplet coding
for the start and end of the AminoAcidRange.

The underlying nucleotide sequence of the AARange can be retrieved
with the labelsubseq() method. This would retrieve the whole DNA
sequence, including possible introns. This is called "DNA_sequence".

To fetch the nucleotide sequence of the Transcript, without introns,
the labelsubseq() of the attached Transcript (the Transcript the
Translation comes from) has to be accessed. This is called
"cDNA_sequence".

Here are the operations to retrieve these latter two kinds of
sequences:

   $startlabel=$AARange->start;
   $endtripletlabel=$AARange->end;
   $endlabel=$AARange->{'seq'}->label(3,$endtripletlabel,$AARange->strand);

   $dnaseq=$AARange->labelsubseq($startlabel,undef,$endlabel));

   $cdnaseq=$AARange->get_Transcript->labelsubseq($startlabel,undef,$endlabel);

To simplify, these operations have been included in two additional
methods: dna_seq() and cdna_seq().

These would return the whole sequence, as in the examples above.  But
the above general scheme can be used by specifying different labels,
to retrieve hypothetical subsequences of interest.

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::AARange;

use strict;
use base qw(Bio::LiveSeq::SeqI);

=head2 new

  Title   : new
  Usage   : $aarange = Bio::LiveSeq::AARange->new(-translation => $obj_ref,
                                               -start => $beginaa,
                                               -end => $endaa,
                                               -name => "ABCD",
                                               -description => "DCBA",
                                               -translength => $length);

  Function: generates a new AminoAcidRange LiveSeq object
  Returns : reference to a new object of class AARange
  Errorcode -1
  Args    : two positions in AminoAcid coordinate numbering
            an object reference specifying to which translation the aminoacid
            ranges refer to
            a name and a description (optional)
            an optional "translength" argument: this can be given when
            a lot of AARanges are to be created at the same time for the same
            Translation object, calculating it with $translation->length
            This would increase the speed, avoiding the new() function to
            calculate everytime the same length again and again for every obj.

=cut

sub new {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my ($obj,%range);

  $obj = \%range;
  $obj = bless $obj, $class;
  my $self=$obj;

  my ($translation,$start,$end,$name,$description,$translength)=($args{-translation},$args{-start},$args{-end},$args{-name},$args{-description},$args{-translength});

  unless (($translation)&&(ref($translation) eq "Bio::LiveSeq::Translation")) {
    $self->warn("No -translation or wrong type given");
    return (-1);
  }
  unless ($translength) { # if it's not given, fetch it
    $translength=$translation->length;
  }
  my $seq=$translation->{'seq'};

  if (($start < 1)&&($start > $translength)) {
    $self->warn("$class not initialised because start aminoacid position not valid");
    return (-1);
  }
  if (($end < 1)&&($end > $translength)) {
    $self->warn("$class not initialised because end aminoacid position not valid");
    return (-1);
  }
  if ($start > $end) {
    $self->warn("$class not initialised because start position > end position!");
    return (-1);
  }

  my ($starttripletlabel,$endtripletlabel);
  if ($start == $end) { # trick to increase speed
    $starttripletlabel=$endtripletlabel=$translation->label($start);
  } else {
    ($starttripletlabel,$endtripletlabel)=($translation->label($start),$translation->label($end));
  }
  unless (($starttripletlabel > 0)&&($endtripletlabel > 0)) {
    $self->warn("$class not initialised because of problems in retrieving start or end label!");
    return (-1);
  }

  # unsure if needed:
  #my $endlabel=$seq->label(3,$endtripletlabel); # to get the real end
  #unless ($endlabel > 0) {
    #carp "$class not initialised because of problems retrieving the last nucleotide of the triplet coding for the end aminoacid";
    #return (-1);
  #}
  $self->{'seq'}=$seq;
  $self->{'start'}=$starttripletlabel;
  $self->{'end'}=$endtripletlabel;
  $self->{'strand'}=$translation->strand;
  $self->{'translation'}=$translation;
  $self->{'name'}=$name;
  $self->{'description'}=$description;
  $self->{'alphabet'}="protein";

  return $obj;
}

sub coordinate_start {
  my $self=shift;
  $self->warn("Cannot perform this operation in an AminoAcidRange object!");
  return (-1);
}

sub all_labels {
  my $self=shift;
  $self->warn("Cannot perform this operation in an AminoAcidRange object!");
  return (-1);
}

sub valid {
  my $self=shift;
  $self->warn("Cannot perform this operation in an AminoAcidRange object!");
  return (-1);
}

=head2 get_Transcript

  Title   : valid
  Usage   : $transcript = $obj->get_Transcript()
  Function: retrieves the reference to the object of class Transcript (if any)
            attached to a LiveSeq object
  Returns : object reference
  Args    : none

=cut

sub get_Transcript {
  my $self=shift;
  return ($self->get_Translation->get_Transcript);
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
  return ($self->{'translation'});
}

sub change {
  my $self=shift;
  $self->warn("Cannot change an AminoAcidRange object!");
  return (-1);
}
sub positionchange {
  my $self=shift;
  $self->warn("Cannot change an AminoAcidRange object!");
  return (-1);
}
sub labelchange {
  my $self=shift;
  $self->warn("Cannot change an AminoAcidRange object!");
  return (-1);
}

sub subseq {
  my ($self,$pos1,$pos2,$length) = @_;
  if (defined ($length)) {
    if ($length < 1) {
      $self->warn("No sense asking for a subseq of length < 1");
      return (-1);
    }
  }
  unless (defined ($pos1)) {
    $pos1=1;
  } elsif ($pos1 < 1) { # if position out of boundaries
    $self->warn("Starting position for AARange cannot be < 1!"); return (-1);
    if ((defined ($pos2))&&($pos1>$pos2)) {
      $self->warn("1st position($pos1) cannot be > 2nd position($pos2)!"); return (-1);
    }
  }
  my $seq=$self->seq;
  my $objlength=length($seq);
  unless (defined ($length)) {
    $length=$objlength-$pos1+1;
  }
  if (defined ($pos2)) {
    if ($pos2 > $objlength) { # if position out of boundaries
      $self->warn("Ending position for AARange cannot be > length of AARange!"); return (-1);
    }
    $length=$pos2-$pos1+1;
    if ((defined ($pos1))&&($pos1>$pos2)) {
      $self->warn("1st position($pos1) cannot be > 2nd position($pos2)!"); return (-1);
    }
  }
  my $str=substr($seq,$pos1-1,$length);
  if (length($str) < $length) {
    $self->warn("Attention, cannot return the length requested for subseq",1);
  }
  return $str;
}

sub seq {
  my $self=shift;
  my ($aa_start,$aa_end)=($self->aa_start,$self->aa_end);
  unless (($aa_start)&&($aa_end)) { # they must both exist
    $self->warn("Not able to find start or end of the AminoAcid Range");
    return (0);
  }
  my $translseq=$self->get_Translation->seq;
  return substr($translseq,$aa_start-1,$aa_end-$aa_start+1);
  # Note: it will return "undef" if the translation stops before the start
  # of the aarange (because of upstream nonsense mutation creating STOP).
  # For the same reason it would return uncomplete (up to the STOP) string
  # if the stop happens in between aarange's start and stop
}

sub length {
  my $self=shift;
  my $seq=$self->seq;
  my $length=length($seq);
  return $length;
}

sub label {
  my ($self,$position)=@_;
  my $translation=$self->get_Translation;
  my $origstart=$translation->coordinate_start; # preserve it
  $translation->coordinate_start($self->start); # change it
  my $label=$translation->label($position);
  $translation->coordinate_start($origstart); # restore it
  return ($label);
}

sub position {
  my ($self,$label)=@_;
  my $translation=$self->get_Translation;
  my $origstart=$translation->coordinate_start; # preserve it
  $translation->coordinate_start($self->start); # change it
  my $position=$translation->position($label);
  $translation->coordinate_start($origstart); # restore it
  return ($position);
}

=head2 aa_start

  Title   : aa_start
  Usage   : $end = $aarange->aa_start()
  Returns : integer (position, according to Translation coordinate system) of
            the start of an AminoAcidRange object
  Args    : none

=cut

sub aa_start {
  my $self=shift;
  my $aastart=$self->get_Translation->position($self->{'start'});
}

=head2 aa_end

  Title   : aa_end
  Usage   : $end = $aarange->aa_end()
  Returns : integer (position, according to Translation coordinate system) of
            the end of an AminoAcidRange object
  Args    : none

=cut

sub aa_end {
  my $self=shift;
  my $aastart=$self->get_Translation->position($self->{'end'});
}

=head2 dna_seq

  Title   : dna_seq
  Usage   : $end = $aarange->dna_seq()
  Returns : the sequence at DNA level of the entire AminoAcidRange
            this would include introns (if present)
  Args    : none

=cut

sub dna_seq {
  my $self=shift;
  my $startlabel=$self->start;
  my $endtripletlabel=$self->end;
  my $endlabel=$self->{'seq'}->label(3,$endtripletlabel,$self->strand);
  return ($self->labelsubseq($startlabel,undef,$endlabel));
}

=head2 cdna_seq

  Title   : cdna_seq
  Usage   : $end = $aarange->cdna_seq()
  Returns : the sequence at cDNA level of the entire AminoAcidRange
            i.e. this is the part of the Transcript that codes for the
            AminoAcidRange. It would be composed just of exonic DNA.
  Args    : none

=cut

sub cdna_seq {
  my $self=shift;
  my $startlabel=$self->start;
  my $endtripletlabel=$self->end;
  my $endlabel=$self->{'seq'}->label(3,$endtripletlabel,$self->strand);
  return ($self->get_Transcript->labelsubseq($startlabel,undef,$endlabel));
}

# this checks if the attached Transcript has a Gene object attached
sub gene {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'gene'} = $value;
  }
  unless (exists $self->{'gene'}) {
    unless (exists $self->get_Transcript->{'gene'}) {
      return (0);
    } else {
      return ($self->get_Transcript->{'gene'});
    }
  } else {
    return $self->{'gene'};
  }
}

1;
