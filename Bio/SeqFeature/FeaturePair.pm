
#
# BioPerl module for Bio::SeqFeature::FeaturePair
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::FeaturePair - hold pair feature information e.g. blast hits

=head1 SYNOPSIS

    my $feat  = new Bio::SeqFeature::FeaturePair(-feature1 => $f1,
						 -feature2 => $f2,
					      );

    # Bio::SeqFeatureI methods can be used

    my $start = $feat->start;
    my $end   = $feat->end;

    # Bio::FeaturePair methods can be used
    my $hstart = $feat->hstart;
    my $hend   = $feat->hend;

   my $feature1 = $feat->feature1;  # returns feature1 object

=head1 DESCRIPTION

A sequence feature object where the feature is itself a feature on another 
sequence - e.g. a blast hit where residues 1-40 of a  protein sequence SW:HBA_HUMAN  
has hit to bases 100 - 220 on a genomic sequence HS120G22.  The genomic sequence 
coordinates are used to create one sequence feature $f1 and the protein coordinates
are used to create feature $f2.  A FeaturePair object can then be made

    my $fp = new Bio::SeqFeature::FeaturePair(-feature1 => $f1,   # genomic
					      -feature2 => $f2,   # protein
					      );

This object can be used as a standard Bio::SeqFeatureI in which case

    my $gstart = $fp->start  # returns start coord on feature1 - genomic seq.
    my $gend   = $fp->end    # returns end coord on feature1.

In general standard Bio::SeqFeatureI method calls return information
in feature1.

Data in the feature 2 object are generally obtained using the standard
methods prefixed by h (for hit!)

    my $pstart = $fp->hstart # returns start coord on feature2 = protein seq.
    my $pend   = $fp->hend   # returns end coord on feature2.

If you wish to swap feature1 and feature2 around :

    $feat->invert

so... 

    $feat->start # etc. returns data in $feature2 object


No sub_SeqFeatures or tags can be stored in this object directly.
Any features or tags are expected to be stored in the contained objects
feature1, and feature2.

=head1 CONTACT

Ewan Birney <birney@sanger.ac.uk>

=head1 DEVELOPERS

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::FeaturePair;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::SeqFeatureI;


@ISA = qw(Bio::Root::Object Bio::SeqFeatureI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize;

  my ($feature1,$feature2) = 
      $self->_rearrange([qw(FEATURE1
			    FEATURE2
			    )],@args);

  # Store the features in the object

  $feature1 && $self->feature1($feature1);
  $feature2 && $self->feature2($feature2);

  # set stuff in self from @args
  return $make; # success - we hope!
}


=head2 feature1

 Title   : feature1
 Usage   : $f = $featpair->feature1
           $featpair->feature1($feature)
 Function: Get/set for the query feature
 Returns : Bio::SeqFeatureI
 Args    : none


=cut


sub feature1 {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("Argument [$arg] must be a Bio::SeqFeatureI") unless (ref($arg) ne "" && $arg->isa("Bio::SeqFeatureI"));
	$self->{_feature1} = $arg;
    } 

    return $self->{_feature1};
}

=head2 feature2

 Title   : feature2
 Usage   : $f = $featpair->feature2
           $featpair->feature2($feature)
 Function: Get/set for the hit feature
 Returns : Bio::SeqFeatureI
 Args    : none


=cut

sub feature2 {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->throw("Argument [$arg] must be a Bio::SeqFeatureI") unless (ref($arg) ne "" && $arg->isa("Bio::SeqFeatureI"));
	$self->{_feature2} = $arg;
    } 
    return $self->{_feature2};
}

=head2 start

 Title   : start
 Usage   : $start = $featpair->start
           $featpair->start(20)
 Function: Get/set on the start coordinate of feature1
 Returns : integer
 Args    : none

=cut

sub start {
    my ($self,$value) = @_;
    
    if (defined($value)) {
	return $self->feature1->start($value);
    } else {
	return $self->feature1->start;
    }

}

=head2 end

 Title   : end
 Usage   : $end = $featpair->end
           $featpair->end($end)
 Function: get/set on the end coordinate of feature1
 Returns : integer
 Args    : none


=cut

sub end{
    my ($self,$value) = @_;

    if (defined($value)) {
	return $self->feature1->end($value);
    } else {
	return $self->feature1->end;
    }
}

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub length {
   my ($self) = @_;

   return $self->end - $self->start +1;
}


=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub strand{
    my ($self,$arg) = @_;

    if (defined($arg)) {
	return $self->feature1->strand($arg);
    } else {
	return $self->feature1->strand;
    }
}

=head2 score

 Title   : score
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub score {
    my ($self,$arg) = @_;
  
    if (defined($arg)) {
	return $self->feature1->score($arg);
    } else {
	return $self->feature1->score;
    }
}

=head2 frame

 Title   : frame
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub frame {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	return $self->feature1->frame($arg);
    } else {
	return $self->feature1->frame;
    }
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $ptag = $featpair->primary_tag
 Function: get/set on the primary_tag of feature1
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub primary_tag{
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	return $self->feature1->primary_tag($arg);
    } else {
	return $self->feature1->primary_tag;
    }
}

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none


=cut

sub source_tag{
    my ($self,$arg) = @_;

    if (defined($arg)) {
	return $self->feature1->source_tag($arg);
    } else {
	return $self->feature1->source_tag;
    }
}

=head2 seqname

 Title   : seqname
 Usage   : $obj->seqname($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store 
           the seqname.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seqname
 Args    : newvalue (optional)


=cut

sub seqname{
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	return $self->feature1->seqname($arg);
    } else {
	return $self->feature1->seqname;
    }
}

=head2 hseqname

 Title   : hseqname
 Usage   : $featpair->hseqname($newval)
 Function: Get/set method for the name of
           feature2.
 Returns : value of $feature2->seqname
 Args    : newvalue (optional)


=cut

sub hseqname {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->feature2->seqname($arg);
    }

    return $self->feature2->seqname;
}


=head2 hstart

 Title   : hstart
 Usage   : $start = $featpair->hstart
           $featpair->hstart(20)
 Function: Get/set on the start coordinate of feature2
 Returns : integer
 Args    : none

=cut

sub hstart {
    my ($self,$value) = @_;
    
    if (defined($value)) {
	return $self->feature2->start($value);
    } else {
	return $self->feature2->start;
    }

}

=head2 hend

 Title   : hend
 Usage   : $end = $featpair->hend
           $featpair->hend($end)
 Function: get/set on the end coordinate of feature2
 Returns : integer
 Args    : none


=cut

sub hend{
    my ($self,$value) = @_;

    if (defined($value)) {
	return $self->feature2->end($value);
    } else {
	return $self->feature2->end;
    }
}



=head2 hstrand

 Title   : hstrand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub hstrand{
    my ($self,$arg) = @_;

    if (defined($arg)) {
	return $self->feature2->strand($arg);
    } else {
	return $self->feature2->strand;
    }
}

=head2 hscore

 Title   : hscore
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub hscore {
    my ($self,$arg) = @_;
  
    if (defined($arg)) {
	return $self->feature2->score($arg);
    } else {
	return $self->feature2->score;
    }
}

=head2 hframe

 Title   : hframe
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub hframe {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	return $self->feature2->frame($arg);
    } else { 
	return $self->feature2->frame;
    }
}

=head2 hprimary_tag

 Title   : hprimary_tag
 Usage   : $ptag = $featpair->hprimary_tag
 Function: Get/set on the primary_tag of feature2
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub hprimary_tag{
    my ($self,$arg) = @_;

    if (defined($arg)) {
	return $self->feature2->primary_tag($arg);
    } else {
	return $self->feature2->primary_tag;
    }
}

=head2 hsource_tag

 Title   : hsource_tag
 Usage   : $tag = $feat->hsource_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none


=cut

sub hsource_tag{
    my ($self,$arg) = @_;

    if (defined($arg)) {
	return $self->feature2->source_tag($arg);
    } else {
	return $self->feature2->source_tag;
    }
}

=head2 invert

 Title   : invert
 Usage   : $tag = $feat->invert
 Function: Swaps feature1 and feature2 around
 Returns : Nothing
 Args    : none


=cut

sub invert {
    my ($self) = @_;

    my $tmp = $self->feature1;
    
    $self->feature1($self->feature2);
    $self->feature2($tmp);

}

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : Function just for complying with SeqFeatureI
 Function:
 Example :
 Returns : an empty list
 Args    :


=cut

sub sub_SeqFeature{
   my ($self,@args) = @_;

   return ();
}

=head2 all_tags

 Title   : all_tags
 Usage   : Function just for complying with SeqFeatureI
 Function:
 Example :
 Returns : an empty list
 Args    :


=cut

sub all_tags{
   my ($self,@args) = @_;
   return ();

}








