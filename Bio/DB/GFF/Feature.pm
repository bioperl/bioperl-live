=head1 NAME

Bio::DB::GFF::Feature -- A relative segment identified by a feature type

=head1 SYNOPSIS

See L<Bio::DB::GFF>.

=head1 DESCRIPTION

Bio::DB::GFF::Feature is a stretch of sequence that corresponding to a
single annotation in a GFF database.  It inherits from
Bio::DB::GFF::RelSegment, and so has all the support for relative
addressing of this class and its ancestors.

Bio::DB::GFF::Feature adds new methods to retrieve the annotation's
type, group, and other GFF attributes.  Annotation types are
represented by Bio::DB::GFF::Typename objects, a simple class that has 
two methods called method() and source().  These correspond to the
method and source fields of a GFF file.  

Annotation groups serve the dual purpose of giving the annotation a
human-readable name, and providing higher-order groupings of
subfeatures into features.  The groups returned by this module are
objects of the Bio::DB::GFF::Featname class.

Bio::DB::GFF::Feature inherits from and implements the abstract
methods of Bio::SeqFeatureI, allowing it to interoperate with other
Bioperl modules.

Generally, you will not create or manipulate Bio::DB::GFF::Feature
objects directly, but use those that are returned by the
Bio::DB::GFF::RelSegment->features() method.

=head1 API

The remainder of this document describes the public and private
methods implemented by this module.

=cut

package Bio::DB::GFF::Feature;

use strict;

use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::GFF::RelSegment;
use Bio::DB::GFF::Featname;
use Bio::DB::GFF::Typename;
use Bio::DB::GFF::Homol;
use Bio::SeqFeatureI;
use Bio::Root::RootI;

use vars qw($VERSION @ISA $AUTOLOAD);
@ISA = qw(Bio::DB::GFF::RelSegment Bio::SeqFeatureI Bio::Root::RootI);

$VERSION = '0.20';

*segments = \&sub_SeqFeature;
*name     = \&group;

=head2 new_from_parent

 Title   : new_from_parent
 Usage   : $f = Bio::DB::GFF::Feature->new_from_parent(@args);
 Function: create a new feature object
 Returns : new Bio::DB::GFF::Feature object
 Args    : see below
 Status  : Internal

This method is called by Bio::DB::GFF to create a new feature using
information obtained from the GFF database.  It is one of two similar
constructors.  This one is called when the feature is generated from a
RelSegment object, and should inherit that object's coordinate system.

The 10 arguments are positional:

  $parent       a Bio::DB::GFF::RelSegment object (or descendent)
  $start        start of this feature
  $stop         stop of this feature
  $method       this feature's GFF method
  $source       this feature's GFF source
  $score	this feature's score
  $fstrand      this feature's strand (relative to the source
                      sequence, which has its own strandedness!)
  $phase        this feature's phase
  $group        this feature's group (a Bio::DB::GFF::Featname object)
  $db_id        this feature's internal database ID                    

=cut

# this is called for a feature that is attached to a parent sequence,
# in which case it inherits its coordinate reference system and strandedness
sub new_from_parent {
  my $package   = shift;
  my ($parent,
      $start,$stop,
      $method,$source,$score,
      $fstrand,$phase,
      $group,$db_id) = @_;

  ($start,$stop) = ($stop,$start) if defined($fstrand) and $fstrand eq '-';
  my $class = $group ? $group->class : $parent->class;

  my $self =  bless {
		     %$parent,
		     start  => $start,
		     stop   => $stop,
		     type   => Bio::DB::GFF::Typename->new($method,$source),
		     fstrand => $fstrand,
		     score  => $score,
		     phase  => $phase,
		     group  => $group,
		     db_id  => $db_id,
		     class  => $class,
		    },$package;
  delete $self->{whole};  # a feature is never whole (sorry)
  $self;
}


=head2 new

 Title   : new
 Usage   : $f = Bio::DB::GFF::Feature->new(@args);
 Function: create a new feature object
 Returns : new Bio::DB::GFF::Feature object
 Args    : see below
 Status  : Internal

This method is called by Bio::DB::GFF to create a new feature using
information obtained from the GFF database.  It is one of two similar
constructors.  This one is called when the feature is generated
without reference to a RelSegment object, and should therefore use its
default coordinate system (relative to itself).

The 11 arguments are positional:

  $factory      a Bio::DB::GFF adaptor object (or descendent)
  $srcseq       the source sequence
  $start        start of this feature
  $stop         stop of this feature
  $method       this feature's GFF method
  $source       this feature's GFF source
  $score	this feature's score
  $fstrand      this feature's strand (relative to the source
                      sequence, which has its own strandedness!)
  $phase        this feature's phase
  $group        this feature's group
  $db_id        this feature's internal database ID                    

=cut

# This is called when creating a feature from scratch.  It does not have
# an inherited coordinate system.
sub new {
  my $package = shift;
  my ($factory,
      $srcseq,
      $start,$stop,
      $method,$source,
      $score,$fstrand,$phase,
      $group,
      $db_id
     ) = @_;

  my $self = bless { },$package;
  ($start,$stop) = ($stop,$start) if defined($fstrand) and $fstrand eq '-';

  my $class = $group ? $group->class : 'Sequence';

  @{$self}{qw(factory sourceseq start stop strand class)} =
    ($factory,$srcseq,$start,$stop,$fstrand,$class);

  @{$self}{qw(ref refstart refstrand)} = ($srcseq,1,$fstrand);

  @{$self}{qw(type fstrand score phase group db_id)} =
    (Bio::DB::GFF::Typename->new($method,$source),$fstrand,$score,$phase,$group,$db_id);

  $self;
}

=head2 type

 Title   : type
 Usage   : $type = $f->type([$newtype])
 Function: get or set the feature type
 Returns : a Bio::DB::GFF::Typename object
 Args    : a new Typename object (optional)
 Status  : Public

This method gets or sets the type of the feature.  The type is a
Bio::DB::GFF::Typename object, which encapsulates the feature method
and source.  

The method() and source() methods described next provide shortcuts to
the individual fields of the type.

=cut

sub type   {
  my $self = shift;
  my $d = $self->{type};
  $self->{type} = shift if @_;
  $d;
}

=head2 method

 Title   : method
 Usage   : $method = $f->method([$newmethod])
 Function: get or set the feature method
 Returns : a string
 Args    : a new method (optional)
 Status  : Public

This method gets or sets the feature method.  It is a convenience
feature that delegates the task to the feature's type object.

=cut

sub method {
  my $self = shift;
  my $d = $self->{type}->method;
  $self->{type}->method(shift) if @_;
  $d;
}

=head2 source

 Title   : source
 Usage   : $source = $f->source([$newsource])
 Function: get or set the feature source
 Returns : a string
 Args    : a new source (optional)
 Status  : Public

This method gets or sets the feature source.  It is a convenience
feature that delegates the task to the feature's type object.

=cut

sub source {
  my $self = shift;
  my $d = $self->{type}->source;
  $self->{type}->source(shift) if @_;
  $d;
}

=head2 score

 Title   : score
 Usage   : $score = $f->score([$newscore])
 Function: get or set the feature score
 Returns : a string
 Args    : a new score (optional)
 Status  : Public

This method gets or sets the feature score.

=cut

sub score  {
  my $self = shift;
  my $d    = $self->{score};
  $self->{score} = shift if @_;
  $d;
}

=head2 phase

 Title   : phase
 Usage   : $phase = $f->phase([$phase])
 Function: get or set the feature phase
 Returns : a string
 Args    : a new phase (optional)
 Status  : Public

This method gets or sets the feature phase.

=cut

sub phase  {
  my $self = shift;
  my $d    = $self->{phase};
  $self->{phase} = shift if @_;
  $d;
}

=head2 strand

 Title   : strand
 Usage   : $strand = $f->strand
 Function: get the feature strand
 Returns : +1, 0 -1
 Args    : none
 Status  : Public

Returns the strand of the feature.  Unlike the other methods, the
strand cannot be changed once the object is created (due to coordinate
considerations).

=cut

sub strand {
  my $self = shift;
  return 0 unless $self->{fstrand};
  return 0 unless defined $self->{start};
  return $self->{start} < $self->{stop} ? '+1' : '-1';
}

=head2 group

 Title   : group
 Usage   : $group = $f->group([$new_group])
 Function: get or set the feature group
 Returns : a Bio::DB::GFF::Featname object
 Args    : a new group (optional)
 Status  : Public

This method gets or sets the feature group.  The group is a
Bio::DB::GFF::Featname object, which has an ID and a class.

=cut

sub group  {
  my $self = shift;
  my $d    = $self->{group};
  $self->{group} = shift if @_;
  $d;
}

=head2 info

 Title   : info
 Usage   : $info = $f->info([$new_info])
 Function: get or set the feature group
 Returns : a Bio::DB::GFF::Featname object
 Args    : a new group (optional)
 Status  : Public

This method is an alias for group().  It is provided for AcePerl
compatibility.

=cut

*info   = \&group;

=head2 target

 Title   : target
 Usage   : $target = $f->target([$new_target])
 Function: get or set the feature target
 Returns : a Bio::DB::GFF::Featname object
 Args    : a new group (optional)
 Status  : Public

This method is an alias for group().  It is provided for AcePerl
compatibility.  In the case of "similarity" features, the group
typically contains the ID of the target sequence, hence the alias.

=cut

*target = \&group;

=head2 id

 Title   : id
 Usage   : $id = $f->id
 Function: get the feature ID
 Returns : a database identifier
 Args    : none
 Status  : Public

This method retrieves the database identifier for the feature.  It
cannot be changed.

=cut

sub id     { shift->{db_id}   }

=head2 clone

 Title   : clone
 Usage   : $feature = $f->clone
 Function: make a copy of the feature
 Returns : a new Bio::DB::GFF::Feature object
 Args    : none
 Status  : Public

This method returns a copy of the feature.

=cut

sub clone {
  my $self = shift;
  my $clone = $self->SUPER::clone;

  if (ref(my $t = $clone->type)) {
    my $type = $t->can('clone') ? $t->clone : bless {%$t},ref $t;
    $clone->type($type);
  }

  if (ref(my $g = $clone->group)) {
    my $group = $g->can('clone') ? $g->clone : bless {%$g},ref $g;
    $clone->group($group);
  }

  $clone;
}



=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feat = $feature->sub_SeqFeature([$method])
 Function: get merged subfeatures
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : a feature method (optional)
 Status  : Public

This method acts like sub_SeqFeature, except that it merges
overlapping segments of the same time into contiguous features.  For
those features that contain heterogeneous subfeatures, you can
retrieve a subset of the subfeatures by providing a method name to
filter on.

=cut

sub sub_SeqFeature {
  my $self = shift;
  my $type = shift;
  my $subfeat = $self->{subfeatures} or return;
  $self->sort_features;
  if ($type) {
    my $features = $subfeat->{lc $type} or return;
    return @{$features};
  } else {
    return map {@{$_}} values %{$subfeat};
  }
}

=head2 add_subfeature

 Title   : add_subfeature
 Usage   : $feature->add_subfeature($feature)
 Function: add a subfeature to the feature
 Returns : nothing
 Args    : a Bio::DB::GFF::Feature object
 Status  : Public

This method adds a new subfeature to the object.  It is used
internally by aggregators, but is available for public use as well.

=cut

sub add_subfeature {
  my $self    = shift;
  my $feature = shift;
  my $type = $feature->method;
  my $subfeat = $self->{subfeatures}{lc $type} ||= [];
  push @{$subfeat},$feature;
}

=head2 merged_segments

 Title   : merged_segments
 Usage   : @segs = $feature->merged_segments([$method])
 Function: get merged subfeatures
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : a feature method (optional)
 Status  : Public

This method acts like sub_SeqFeature, except that it merges
overlapping segments of the same time into contiguous features.  For
those features that contain heterogeneous subfeatures, you can
retrieve a subset of the subfeatures by providing a method name to
filter on.

A side-effect of this method is that the features are returned in
sorted order by their start tposition.

=cut

sub merged_segments {
  my $self = shift;
  my $type = shift;
  $type ||= '';    # prevent uninitialized variable warnings

  return @{$self->{merged_segs}{$type}} if exists $self->{merged_segs}{$type};

  my @segs = sort {$a->start <=> $b->start} $self->sub_SeqFeature($type);
  # attempt to merge overlapping segments
  my @merged;
  for my $s (@segs) {
    my $previous = $merged[-1];
    if ($previous && $previous->stop+1 >= $s->start) {
      $previous->{stop} = $s->{stop};
      # fix up the target too
      my $g = $previous->{group};
      if ( ref($g) &&  $g->isa('Bio::DB::GFF::Homol')) {
	my $cg = $s->{group};
	$g->{stop} = $cg->{stop};
      }
    } else {
      my $copy = $s->clone;
      push @merged,$copy;
    }
  }
  $self->{merged_segs}{$type} = \@merged;
  return @merged;
}

=head2 sub_types

 Title   : sub_types
 Usage   : @methods = $feature->sub_types
 Function: get methods of all sub-seqfeatures
 Returns : a list of method names
 Args    : none
 Status  : Public

For those features that contain subfeatures, this method will return a
unique list of method names of those subfeatures, suitable for use
with sub_SeqFeature().

=cut

sub sub_types {
  my $self = shift;
  my $subfeat = $self->{subfeatures} or return;
  return keys %$subfeat;
}

=head2 Autogenerated Methods

 Title   : AUTOLOAD
 Usage   : @subfeat = $feature->Method
 Function: Return subfeatures using autogenerated methods
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : none
 Status  : Public

Any method that begins with an initial capital letter will be passed
to AUTOLOAD and treated as a call to sub_SeqFeature with the method
name used as the method argument.  For instance, this call:

  @exons = $feature->Exon;

is equivalent to this call:

  @exons = $feature->sub_SeqFeature('exon');

=cut

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  my $sub = $AUTOLOAD;
  my $self = $_[0];

  # ignore DESTROY calls
  return if $func_name eq 'DESTROY';

  # fetch subfeatures if func_name has an initial cap
  return $self->sub_SeqFeature($func_name) if $func_name =~ /^[A-Z]/;

  # error message of last resort
  $self->throw(qq(Can\'t locate object method "$func_name" via package "$pack"));
}

=head2 adjust_bounds

 Title   : adjust_bounds
 Usage   : $feature->adjust_bounds
 Function: adjust the bounds of a feature
 Returns : ($start,$stop,$strand)
 Args    : none
 Status  : Public

This method adjusts the boundaries of the feature to enclose all its
subfeatures.  It returns the new start, stop and strand of the
enclosing feature.

=cut

# adjust a feature so that its boundaries are synched with its subparts' boundaries.
# this works recursively, so subfeatures can contain other features
sub adjust_bounds {
  my $self = shift;
  my $g = $self->{group};

  if (my $subfeat = $self->{subfeatures}) {
    for my $list (values %$subfeat) {
      for my $feat (@$list) {

	# fix up our bounds to hold largest subfeature
	my($start,$stop,$strand) = $feat->adjust_bounds;
	$self->{fstrand} = $strand unless defined $self->{fstrand};
	if ($self->strand >= 0) {
	  $self->{start} = $start if !defined($self->{start}) || $start < $self->{start};
	  $self->{stop}  = $stop  if !defined($self->{stop})  || $stop  > $self->{stop};
	} else {
	  $self->{start} = $start if !defined($self->{start}) || $start > $self->{start};
	  $self->{stop}  = $stop  if !defined($self->{stop})  || $stop  < $self->{stop};
	}

	# fix up endpoints of targets too (for homologies only)
	my $h = $feat->group;
	next unless $h && $h->isa('Bio::DB::GFF::Homol');
	next unless $g && $g->isa('Bio::DB::GFF::Homol');
	($start,$stop) = ($h->{start},$h->{stop});
	if ($h->strand >= 0) {
	  $g->{start} = $start if !defined($g->{start}) || $start < $g->{start};
	  $g->{stop}  = $stop  if !defined($g->{stop})  || $stop  > $g->{stop};
	} else {
	  $g->{start} = $start if !defined($g->{start}) || $start > $g->{start};
	  $g->{stop}  = $stop  if !defined($g->{stop})  || $stop  < $g->{stop};
	}
      }
    }
  }

  ($self->{start},$self->{stop},$self->strand);
}

=head2 sort_features

 Title   : sort_features
 Usage   : $feature->sort_features
 Function: sort features
 Returns : nothing
 Args    : none
 Status  : Public

This method sorts subfeatures in ascending order by their start
position.  For reverse strand features, it sorts subfeatures in
descending order.  After this is called sub_SeqFeature will return the
features in order.

This method is called internally by merged_segments().

=cut

# sort features
sub sort_features {
  my $self = shift;
  return if $self->{sorted}++;
  my $strand = $self->strand or return;
  my $subfeat = $self->{subfeatures} or return;
  for my $type (keys %$subfeat) {
    $subfeat->{$type} = [sort {$a->start<=>$b->start} @{$subfeat->{$type}}] if $strand > 0;
    $subfeat->{$type} = [sort {$b->start<=>$a->start} @{$subfeat->{$type}}] if $strand < 0;
  }
}

=head2 asString

 Title   : asString
 Usage   : $string = $feature->asString
 Function: return human-readabled representation of feature
 Returns : a string
 Args    : none
 Status  : Public

This method returns a human-readable representation of the feature and
is called by the overloaded "" operator.

=cut

sub asString {
  my $self = shift;
  my $type = $self->type;
  my $name = $self->group;
  return "$type($name)" if $name;
  return $type;
#  my $type = $self->method;
#  my $id   = $self->group || 'unidentified';
#  return join '/',$id,$type,$self->SUPER::asString;
}

=head1 A Note About Similarities

The current default aggregator for GFF "similarity" features creates a
composite Bio::DB::GFF::Feature object of type "gapped_alignment".
The target() method for the feature as a whole will return a
RelSegment object that is as long as the extremes of the similarity
hit target, but will not necessarily be the same length as the query
sequence.  The length of each "similarity" subfeature will be exactly
the same length as its target().  These subfeatures are essentially
the HSPs of the match.

The following illustrates this:

  @similarities = $segment->feature('similarity:BLASTN');
  $sim          = $similarities[0];

  print $sim->type;        # yields "gapped_similarity:BLASTN"

  $query_length  = $sim->length;
  $target_length = $sim->target->length;  # $query_length != $target_length

  @matches = $sim->Similarity;   # use autogenerated method
  $query1_length  = $matches[0]->length;
  $target1_length = $matches[0]->target->length; # $query1_length == $target1_length

If you merge segments by calling merged_segments(), then the length of
the query sequence segments will no longer necessarily equal the
length of the targets, because the alignment information will have
been lost.  Nevertheless, the targets are adjusted so that the first
and last base pairs of the query match the first and last base pairs
of the target.

=cut

1;

=head1 BUGS

This module is still under development.

=head1 SEE ALSO

L<bioperl>, L<Bio::DB::GFF>, L<Bio::DB::RelSegment>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

