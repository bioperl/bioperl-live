#
# bioperl module for Bio::LiveSeq::Translation
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

Bio::LiveSeq::Translation - Translation class for LiveSeq

=head1 SYNOPSIS

  #documentation needed

=head1 DESCRIPTION

This stores information about aminoacids translations of transcripts.
The implementation is that a Translation object is the translation of
a Transcript object, with different possibilities of manipulation,
different coordinate system and eventually its own ranges (protein domains).

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::Translation;

use strict;
#use Carp qw(croak carp cluck);
use Bio::LiveSeq::SeqI; # uses SeqI, inherits from it
use Bio::PrimarySeq;
use base qw(Bio::LiveSeq::Transcript);


=head2 new

  Title   : new
  Usage   : $protein = Bio::LiveSeq::Translation->new(-transcript => $transcr);

  Function: generates a new Bio::LiveSeq::Translation
  Returns : reference to a new object of class Translation
  Errorcode -1
  Args    : reference to an object of class Transcript

=cut

sub new {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my ($obj,%translation);

  my $transcript=$args{-transcript};

  $obj = \%translation;
  $obj = bless $obj, $class;

  unless ($transcript) {
    $obj->throw("$class not initialised because no -transcript given");
  }
  unless (ref($transcript) eq "Bio::LiveSeq::Transcript") {
    $obj->throw("$class not initialised because no object of class Transcript given");
  }

  #my $startbase = $transcript->start;
  #my $endbase = $transcript->end;
  my $strand = $transcript->strand;
  my $seq = $transcript->{'seq'};

  $obj->{'strand'}=$strand;
  $obj->{'seq'}=$seq;
  $obj->{'transcript'}=$transcript;
  $obj->{'alphabet'}="protein";

  $transcript->{'translation'}=$obj;# set the Translation ref into its Transcript
  return $obj;
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
  return ($self->{'transcript'});
}

# These get redefined here, overriding the SeqI ones

sub change {
  my ($self)=@_;
  $self->warn("Cannot change a Translation object!\nChanges have to be issued at the nucleotide level!");
  return (-1);
}
sub positionchange {
  my ($self)=@_;
  $self->warn("Cannot change a Translation object!\nChanges have to be issued at the nucleotide level!");
  return (-1);
}
sub labelchange {
  my ($self)=@_;
  $self->warn("Cannot change a Translation object!\nChanges have to be issued at the nucleotide level!");
  return (-1);
}

# this just returns the translation of the transcript, without checking for
# stop codons
sub transl_seq {
  my $self=shift;
  my $transcript=$self->get_Transcript;
  my $translation=$transcript->translate(undef, undef, undef, 
					 $self->translation_table)->seq;
  return $translation;
}

# version 1.74 -> now the "*" is printed
sub seq {
  my $self=shift;
  my $proteinseq;
  my $transcript=$self->get_Transcript;
  my $translation=$transcript->translate(undef, undef, undef, 
					 $self->translation_table)->seq;
  my $stop_pos=index($translation,"*");
  if ($stop_pos == -1) { # no stop present, continue downstream
    my $downstreamseq=$transcript->downstream_seq();
    #carp "the downstream is: $downstreamseq"; # debug
    my $cdnaseq=$transcript->seq();
    my $extendedseq = Bio::PrimarySeq->new(-seq => "$cdnaseq$downstreamseq",
					  -alphabet => 'dna'
					  );

    $translation=$extendedseq->translate(undef, undef, undef, 
					 $self->translation_table)->seq;
    #carp "the new translation is: $translation"; # debug
    $stop_pos=index($translation,"*");
    if ($stop_pos == -1) { # still no stop present, return warning
      $self->warn("Warning: no stop codon found in the retrieved sequence downstream of Transcript ",1);
      undef $stop_pos;
      $proteinseq=$translation;
    } else {
      $proteinseq=substr($translation,0,$stop_pos+1);
      #carp "the new stopped translation is: $proteinseq, because the stop is at position $stop_pos"; # debug
    }
  } else {
    $proteinseq=substr($translation,0,$stop_pos+1);
  }
  return $proteinseq;
}

sub length {
  my $self=shift;
  my $seq=$self->seq;
  my $length=length($seq);
  return $length;
}

sub all_labels {
  my $self=shift;
  return $self->get_Transcript->all_labels;
}

# counts in triplet. Only a label matching the beginning of a triplet coding
# for an aminoacid is considered valid when setting coordinate_start
# (i.e. only in frame!)
sub valid {
  my ($self,$label)=@_;
  my $i;
  my @labels=$self->get_Transcript->all_labels;
  my $length=$#labels;
  while ($i <= $length) {
    if ($label == $labels[$i]) {
      return (1); # found
    }
    $i=$i+3;
  }
  return (0); # not found
}

# returns the label to the first nucleotide of the triplet coding for $position aminoacid
sub label {
  my ($self,$position)=@_;
  my $firstlabel=$self->coordinate_start; # this is in_frame checked
  if ($position > 0) {
    $position=$position*3-2;
  } else { # if position = 0 this will be caught by Transcript, error thrown
    $position=$position*3;
  }
  return $self->get_Transcript->label($position,$firstlabel);
  # check for coord_start different
}

# returns position (aminoacids numbering) of a particular label
# used to return 0 for not in frame labels
# now returns the position anyway (after version 1.66)
sub position {
  my ($self,$label)=@_;
  my $firstlabel=$self->coordinate_start; # this is in_frame checked
  my $position=$self->get_Transcript->position($label,$firstlabel);
  use integer;
  my $modulus=$position % 3;
  if ($position == 0) {
    return (0);
  } elsif ($position > 0) {
    if ($modulus != 1) {
      $self->warn("Attention! Label $label is not in frame ". 
		  "(1st position of triplet) with protein",1) if $self->verbose > 0; # ignorable
      if ($modulus == 2) {
	return ($position / 3 + 1);
      } else { # i.e. modulus == 0
	return ($position / 3);
      }
    }
    return ($position / 3 + 1);
  } else { # pos < 0
    if ($modulus != 0) {
      $self->warn("Attention! Label $label is not in frame ".
		  "(1st position of triplet) with protein",1) if $self->verbose > 0; # ignorable
      return ($position / 3 - 1); # ok for both other positions
    }
    return ($position / 3);
  }
  $self->throw( "WEIRD: execution shouldn't have reached here");
  return (0); # this should never happen, but just in case
}

# note: it inherits subseq and labelsubseq from Transcript!

sub start {
  my $self=shift;
  return ($self->{'transcript'}->start);
}

sub end {
  my $self=shift;
  return ($self->{'transcript'}->end);
}

=head2 aa_ranges

  Title   : aa_ranges
  Usage   : @proteinfeatures = $translation->aa_ranges()
  Function: to retrieve all the LiveSeq AARange objects attached to a
            Translation, usually created out of a SwissProt database entry
            crossreferenced from an EMBL CDS feature.
  Returns : an array
  Args    : none

=cut

# returns an array of obj_ref of AARange objects attached to the Translation
sub aa_ranges {
  my $self=shift;
  return ($self->{'aa_ranges'});
}

sub translation_table {
  my $self=shift;
  $self->get_Transcript->translation_table(@_);
}

# returns all aminoacids "affected" i.e. all aminoacids coded by any codon
# "touched" by the range selected between the labels, even if only partially.

# it's not optimized for performance but it's useful

sub labelsubseq {
  my ($self,$start,$length,$end)=@_;
  my ($pos1,$pos2);
  my $transcript=$self->get_Transcript;
  if ($start) {
    unless ($transcript->valid($start)) {
      $self->warn("Start label not valid"); return (-1);
    }
    $pos1=$self->position($start);
  }
  if ($end) {
    if ($end == $start) {
      $length=1;
    } else {
      unless ($transcript->valid($end)) {
	$self->warn("End label not valid"); return (-1);
      }
      unless ($transcript->follows($start,$end) == 1) {
	$self->warn("End label does not follow Start label!"); return (-1);
      }
      $pos2=$self->position($end);
      $length=$pos2-$pos1+1;
    }
  }
  my $sequence=$self->seq;
  return (substr($sequence,$pos1-1,$length));
}

# return the offset in aminoacids from LiveSeq protein sequence and SwissProt
# sequence (usually as a result of an INIT_MET or a gap)
sub offset {
  my $self=shift;
  return ($self->{'offset'});
}

1;
