# $Id$
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

A sequence feature object where the feature is itself a feature on
another sequence - e.g. a blast hit where residues 1-40 of a protein
sequence SW:HBA_HUMAN has hit to bases 100 - 220 on a genomic sequence
HS120G22.  The genomic sequence coordinates are used to create one
sequence feature $f1 and the protein coordinates are used to create
feature $f2.  A FeaturePair object can then be made

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


No sub_SeqFeatures or tags can be stored in this object directly.  Any
features or tags are expected to be stored in the contained objects
feature1, and feature2.

=head1 CONTACT

Ewan Birney E<lt>birney@sanger.ac.ukE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::FeaturePair;
use vars qw(@ISA);
use strict;

use Bio::SeqFeatureI;
use Bio::SeqFeature::Generic;

@ISA = qw(Bio::SeqFeature::Generic);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($feature1,$feature2) = 
	$self->_rearrange([qw(FEATURE1
			      FEATURE2
			      )],@args);
    
    # Store the features in the object
    $feature1 && $self->feature1($feature1);
    $feature2 && $self->feature2($feature2);
    return $self;
}

=head2 feature1

 Title   : feature1
 Usage   : $f = $featpair->feature1
           $featpair->feature1($feature)
 Function: Get/set for the query feature
 Returns : Bio::SeqFeatureI
 Args    : Bio::SeqFeatureI


=cut

sub feature1 {
    my ($self,$arg) = @_;    
    if ( defined($arg) || !defined $self->{'feature1'} ) {
	$arg = new Bio::SeqFeature::Generic() unless( defined $arg);
	$self->throw("Argument [$arg] must be a Bio::SeqFeatureI") 
	    unless (ref($arg) && $arg->isa("Bio::SeqFeatureI"));
	$self->{'feature1'} = $arg;
    }
    return $self->{'feature1'};
}

=head2 feature2

 Title   : feature2
 Usage   : $f = $featpair->feature2
           $featpair->feature2($feature)
 Function: Get/set for the hit feature
 Returns : Bio::SeqFeatureI
 Args    : Bio::SeqFeatureI


=cut

sub feature2 {
    my ($self,$arg) = @_;

    if ( defined($arg) || ! defined $self->{'feature2'}) {
	$arg = new Bio::SeqFeature::Generic unless( defined $arg);
	$self->throw("Argument [$arg] must be a Bio::SeqFeatureI") 
	    unless (ref($arg) && $arg->isa("Bio::SeqFeatureI"));
	$self->{'feature2'} = $arg;
    }
    return $self->{'feature2'};
}

# Internal overridable getter/setter for the actual stored value of
# seq_id.  Delegates to the same method in feature1.
sub _seq_id {
  my $self = shift;
  $self->feature1()->_seq_id( @_ );
} # _seq_id(..)

# Internal overridable getter/setter for the actual stored value of
# start.  Delegates to the same method in feature1.
sub _start {
  my $self = shift;
  $self->feature1()->_start( @_ );
} # _start(..)

# Internal overridable getter/setter for the actual stored value of
# end.  Delegates to the same method in feature1.
sub _end {
  my $self = shift;
  $self->feature1()->_end( @_ );
} # _end(..)

# Internal overridable getter/setter for the actual stored value of
# strand.  Delegates to the same method in feature1.
sub _strand {
  my $self = shift;
  $self->feature1()->_strand( @_ );
} # _strand(..)

=head2 location

 Title   : location
 Usage   : $location = $featpair->location
           $featpair->location($location)
 Function: Get/set location object (using feature1)
 Returns : Bio::LocationI object
 Args    : [optional] LocationI to store

=cut

sub location {
    my ($self,$value) = @_;    
    return $self->feature1->location($value);
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
    return $self->feature1->score($arg);    
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
    return $self->feature1->frame($arg);    
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
    return $self->feature1->primary_tag($arg);    
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
    return $self->feature1->source_tag($arg);    
}

=head2 seqname

 Title   : seqname
 Usage   : $obj->seq_id($newval)
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
    return $self->feature1->seq_id($arg);    
}

=head2 hseqname

 Title   : hseqname
 Usage   : $featpair->hseqname($newval)
 Function: Get/set method for the name of
           feature2.
 Returns : value of $feature2->seq_id
 Args    : newvalue (optional)


=cut

sub hseqname {
    my ($self,$arg) = @_;
    return $self->feature2->seq_id($arg);
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
    return $self->feature2->start($value);    
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
    return $self->feature2->end($value);    
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
    return $self->feature2->strand($arg);
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
    return $self->feature2->score($arg);    
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
    return $self->feature2->frame($arg);    
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
    return $self->feature2->primary_tag($arg);    
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
    return $self->feature2->source_tag($arg);
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
    return undef;
}

1;
