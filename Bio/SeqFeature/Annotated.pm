package Bio::SeqFeature::Annotated;

use strict;
use base qw(Bio::Root::Root Bio::SeqFeatureI Bio::AnnotatableI Bio::FeatureHolderI);

use Bio::Root::Root;
use Bio::Annotation::Collection;
use Bio::Location::Simple;
use Bio::Tools::GFF;

sub new {
    my ( $caller, @args) = @_;
    my ($self) = $caller->SUPER::new(@args); 

    $self->_initialize(@args);

    return $self;
}

sub _initialize {
  my ($self,@args) = @_;
  my ($start, $end, $strand, $frame,$score,
      $seq_id, $annot, $location,$display_name) =
        $self->_rearrange([qw(START
                              END
                              STRAND
                              FRAME
                              SCORE
                              SEQ_ID
                              ANNOTATION
                              LOCATION
                              DISPLAY_NAME
                             )], @args);

  defined $start        && $self->start($start);
  defined $end          && $self->end($end);
  defined $strand       && $self->strand($strand);
  defined $frame        && $self->frame($frame);
  defined $display_name && $self->display_name($display_name);
  defined $score        && $self->score($score);
  $location             && $self->location($location);
  $annot                && $self->annotation($annot);
  $seq_id               && $self->seq_id($seq_id);

}

=head2 attach_seq

 Title   : attach_seq
 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000
 Example :
 Returns : TRUE on success
 Args    : a Bio::PrimarySeqI compliant object


=cut

sub attach_seq {
   my ($self, $seq) = @_;

   if ( ! ($seq && ref($seq) && $seq->isa("Bio::PrimarySeqI")) ) {
       $self->throw("Must attach Bio::PrimarySeqI objects to SeqFeatures");
   }

   $self->{'_gsf_seq'} = $seq;

   # attach to sub features if they want it
   foreach ( $self->sub_SeqFeature() ) {
       $_->attach_seq($seq);
   }
   return 1;
}

=head2 seq

 Title   : seq
 Usage   : $tseq = $sf->seq()
 Function: returns the truncated sequence (if there) for this
 Example :
 Returns : sub seq (a Bio::PrimarySeqI compliant object) on attached sequence
           bounded by start & end, or undef if there is no sequence attached
 Args    : none


=cut

sub seq {
   my ($self, $arg) = @_;

   if ( defined $arg ) {
       $self->throw("Calling SeqFeature::Generic->seq with an argument. You probably want attach_seq");
   }

   if ( ! exists $self->{'_gsf_seq'} ) {
       return undef;
   }

   # assumming our seq object is sensible, it should not have to yank
   # the entire sequence out here.

   my $seq = $self->{'_gsf_seq'}->trunc($self->start(), $self->end());


   if ( defined $self->strand &&
	$self->strand == -1 ) {

       # ok. this does not work well (?)
       #print STDERR "Before revcom", $seq->str, "\n";
       $seq = $seq->revcom;
       #print STDERR "After  revcom", $seq->str, "\n";
   }

   return $seq;
}

=head2 entire_seq

 Title   : entire_seq
 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns : a Bio::PrimarySeqI compliant object, or undef if there is no
           sequence attached
 Args    :


=cut

sub entire_seq {
   return shift->{'_gsf_seq'};
}


=head2 seq_id

 Title   : seq_id
 Usage   : $obj->seq_id($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store
           the ID (e.g., display_id) of the sequence.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seq_id
 Args    : newvalue (optional)


=cut

sub seq_id {
    my $obj = shift;
    return $obj->{'_gsf_seq_id'} = shift if @_;
    return $obj->{'_gsf_seq_id'};
}

=head2 display_name

 Title   : display_name
 Usage   : $featname = $obj->display_name
 Function: Implements the display_name() method, which is a human-readable
           name for the feature. 
 Returns : value of display_name (a string)
 Args    : Optionally, on set the new value or undef 

=cut

sub display_name{
    my $self = shift;
    return $self->{'display_name'} = shift if @_;
    return $self->{'display_name'};
}

=head1 Methods for implementing Bio::AnnotatableI

=cut

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($annot_obj)
 Function: Get/set the annotation collection object for annotating this
           feature.

 Example : 
 Returns : A Bio::AnnotationCollectionI object
 Args    : newvalue (optional)


=cut

sub annotation {
    my ($obj,$value) = @_;

    # we are smart if someone references the object and there hasn't been
    # one set yet
    if(defined $value || ! defined $obj->{'annotation'} ) {
        $value = new Bio::Annotation::Collection unless ( defined $value );
        $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};
}

=head1 Methods to implement Bio::FeatureHolderI

This includes methods for retrieving, adding, and removing
features. Since this is already a feature, features held by this
feature holder are essentially sub-features.

=cut

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : @feats = $feat->get_SeqFeatures();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none


=cut

sub get_SeqFeatures {
    return @{ shift->{'_gsf_sub_array'} || []};
}

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $feat->add_SeqFeature($subfeat);
           $feat->add_SeqFeature($subfeat,'EXPAND')
 Function: adds a SeqFeature into the subSeqFeature array.
           with no 'EXPAND' qualifer, subfeat will be tested
           as to whether it lies inside the parent, and throw
           an exception if not.

           If EXPAND is used, the parent's start/end/strand will
           be adjusted so that it grows to accommodate the new
           subFeature
 Returns : nothing
 Args    : An object which has the SeqFeatureI interface


=cut

#'
sub add_SeqFeature{
    my ($self,$feat,$expand) = @_;
    unless( defined $feat ) {
	$self->warn("Called add_SeqFeature with no feature, ignoring");
	return;
    }
    if ( !$feat->isa('Bio::SeqFeatureI') ) {
        $self->warn("$feat does not implement Bio::SeqFeatureI. Will add it anyway, but beware...");
    }

    if($expand && ($expand eq 'EXPAND')) {
        $self->_expand_region($feat);
    } else {
        if ( !$self->contains($feat) ) {
            $self->throw("$feat is not contained within parent feature, and expansion is not valid");
        }
    }

    $self->{'_gsf_sub_array'} = [] unless exists($self->{'_gsf_sub_array'});
    push(@{$self->{'_gsf_sub_array'}},$feat);

}

=head2 remove_SeqFeatures

 Title   : remove_SeqFeatures
 Usage   : $sf->remove_SeqFeatures
 Function: Removes all sub SeqFeatures

           If you want to remove only a subset, remove that subset from the
           returned array, and add back the rest.

 Example :
 Returns : The array of Bio::SeqFeatureI implementing sub-features that was
           deleted from this feature.
 Args    : none


=cut

sub remove_SeqFeatures {
   my ($self) = @_;

   my @subfeats = @{$self->{'_gsf_sub_array'} || []};
   $self->{'_gsf_sub_array'} = []; # zap the array implicitly.
   return @subfeats;
}

=head2 _expand_region

 Title   : _expand_region
 Usage   : $self->_expand_region($feature);
 Function: Expand the total region covered by this feature to
           accomodate for the given feature.

           May be called whenever any kind of subfeature is added to this
           feature. add_sub_SeqFeature() already does this.
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.


=cut

sub _expand_region {
    my ($self, $feat) = @_;
    if(! $feat->isa('Bio::SeqFeatureI')) {
        $self->warn("$feat does not implement Bio::SeqFeatureI");
    }
    # if this doesn't have start/end set - forget it!
    if((! defined($self->start())) && (! defined $self->end())) {
        $self->start($feat->start());
        $self->end($feat->end());
#        $self->strand($feat->strand) unless defined($self->strand());
        $self->strand($feat->strand) unless $self->strand();
    } else {
        my $range = $self->union($feat);
        $self->start($range->start);
        $self->end($range->end);
        $self->strand($range->strand);
    }
}

=head2 source

 Title   : source
 Usage   : $obj->source($newval)
 Function: holds a string corresponding to the source of the feature.
           this method may be moved to a Bio::Annotation::SimpleValue in the future.
 Example : 
 Returns : value of source (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub source {
  my $self = shift;

  return $self->{'source'} = shift if @_;
  return $self->{'source'};
}


=head2 start

 Title   : start
 Usage   : $start = $feat->start
           $feat->start(20)
 Function: Get/set on the start coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub start {
  my ($self,$value) = @_;
  return $self->location->start($value);
}

=head2 end

 Title   : end
 Usage   : $end = $feat->end
           $feat->end($end)
 Function: get/set on the end coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub end {
  my ($self,$value) = @_;
  return $self->location->end($value);
}

=head2 length

 Title   : length
 Usage   : my $len = $feature->length
 Function: Get the feature length computed as 
           $feat->end - $feat->start + 1
 Returns : integer
 Args    : none


=cut

sub length {
  my $self = shift;
  return $self->end() - $self->start() + 1;
}

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub strand {
  my $self = shift;
  return $self->location->strand(@_);
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
  my $self = shift;

  if (@_) {
    my $value = shift;
    if ( defined $value && $value && $value !~ /^[+-]?\d+\.?\d*(e-\d+)?/ and $value != 0) {
	  $self->throw(-class=>'Bio::Root::BadParameter',
                   -text=>"'$value' is not a valid score",
                   -value=>$value);
    }
    return $self->{'_gsf_score'} = $value;
  }
  return $self->{'_gsf_score'};
}

=head2 frame

 Title   : frame
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2, '.'
 Args    : none if get, the new value if set


=cut

sub frame {
  my $self = shift;

  if ( @_ ) {
    my $value = shift;
    if ( defined $value && 
         $value !~ /^[0-2.]$/ ) {
	  $self->throw("'$value' is not a valid frame");
    }
    if( defined $value && $value eq '.' ) { $value = '.' } 
    return $self->{'_gsf_frame'} = $value;
  }
  return $self->{'_gsf_frame'};
}

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location 
           of feature on sequence or parent feature  
 Returns : Bio::LocationI object
 Args    : [optional] Bio::LocationI object to set the value to.


=cut

sub location {
  my($self, $value ) = @_;

  if (defined($value)) {
    unless (ref($value) and $value->isa('Bio::LocationI')) {
      $self->throw("object $value pretends to be a location but ".
                   "does not implement Bio::LocationI");
    }
    $self->{'_location'} = $value;
  }
  elsif (! $self->{'_location'}) {
    # guarantees a real location object is returned every time
    $self->{'_location'} = Bio::Location::Simple->new();
  }
  return $self->{'_location'};
}

sub _no_tags {
  my $self = shift;
  $self->throw("tag methods are deprecated. use Bio::Annotation::Collection");
  #$self->throw_not_implemented();
}

sub primary_tag {
  return shift->_no_tags();
}

sub source_tag {
  return shift->_no_tags();
}

sub has_tag {
  return shift->_no_tags();
}

sub add_tag_value {
  return shift->_no_tags();
}

sub get_tag_values {
  return shift->_no_tags();
}

sub get_all_tags {
  return shift->_no_tags();
}

sub remove_tag {
  return shift->_no_tags();
}

1;
