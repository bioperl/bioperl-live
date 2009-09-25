=head1 NAME

Bio::DB::GFF::Feature -- A relative segment identified by a feature type

=head1 SYNOPSIS

See L<Bio::DB::GFF>.

=head1 DESCRIPTION

Bio::DB::GFF::Feature is a stretch of sequence that corresponding to a
single annotation in a GFF database.  It inherits from
Bio::DB::GFF::RelSegment, and so has all the support for relative
addressing of this class and its ancestors.  It also inherits from
Bio::SeqFeatureI and so has the familiar start(), stop(),
primary_tag() and location() methods (it implements Bio::LocationI
too, if needed).

Bio::DB::GFF::Feature adds new methods to retrieve the annotation
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
Bio::DB::GFF::RelSegment-E<gt>features() method.

=head2 Important note about start() vs end()

If features are derived from segments that use relative addressing
(which is the default), then start() will be less than end() if the
feature is on the opposite strand from the reference sequence.  This
breaks Bio::SeqI compliance, but is necessary to avoid having the real
genomic locations designated by start() and end() swap places when
changing reference points.

To avoid this behavior, call $segment-E<gt>absolute(1) before fetching
features from it.  This will force everything into absolute
coordinates.

For example:

 my $segment = $db->segment('CHROMOSOME_I');
 $segment->absolute(1);
 my @features = $segment->features('transcript');

=head1 API

The remainder of this document describes the public and private
methods implemented by this module.

=cut

package Bio::DB::GFF::Feature;

use strict;

use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::GFF::Featname;
use Bio::DB::GFF::Typename;
use Bio::DB::GFF::Homol;
use Bio::LocationI;
use Data::Dumper;

use vars qw($AUTOLOAD);
use base qw(Bio::DB::GFF::RelSegment Bio::SeqFeatureI Bio::Root::Root);

#' 

*segments = *get_SeqFeatures = \&sub_SeqFeature;

my %CONSTANT_TAGS = (method=>1, source=>1, score=>1, phase=>1, notes=>1, id=>1, group=>1);

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
RelSegment object, and should inherit the coordinate system of that 
object.

The 13 arguments are positional (sorry):

  $parent       a Bio::DB::GFF::RelSegment object (or descendent)
  $start        start of this feature
  $stop         stop of this feature
  $method       this feature's GFF method
  $source       this feature's GFF source
  $score	       this feature's score
  $fstrand      this feature's strand (relative to the source
                      sequence, which has its own strandedness!)
  $phase        this feature's phase
  $group        this feature's group (a Bio::DB::GFF::Featname object)
  $db_id        this feature's internal database ID
  $group_id     this feature's internal group database ID
  $tstart       this feature's target start
  $tstop        this feature's target stop

tstart and tstop are not used for anything at the moment, since the
information is embedded in the group object.

=cut

# this is called for a feature that is attached to a parent sequence,
# in which case it inherits its coordinate reference system and strandedness
sub new_from_parent {
  my $package   = shift;
  my ($parent,
      $start,$stop,
      $method,$source,$score,
      $fstrand,$phase,
      $group,$db_id,$group_id,
      $tstart,$tstop) = @_;

  ($start,$stop) = ($stop,$start) if defined($fstrand) and $fstrand eq '-';
  my $class = $group ? $group->class : $parent->class;

  my $self =  bless {
		     factory   => $parent->{factory},
		     sourceseq => $parent->{sourceseq},
		     strand    => $parent->{strand},
		     ref       => $parent->{ref},
		     refstart  => $parent->{refstart},
		     refstrand => $parent->{refstrand},
		     absolute  => $parent->{absolute},
		     start     => $start,
		     stop      => $stop,
		     type      => Bio::DB::GFF::Typename->new($method,$source),
		     fstrand   => $fstrand,
		     score     => $score,
		     phase     => $phase,
		     group     => $group,
		     db_id     => $db_id,
		     group_id  => $group_id,
		     class     => $class,
		    },$package;
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
  $score	       this feature's score
  $fstrand      this feature's strand (relative to the source
                      sequence, which has its own strandedness!)
  $phase        this feature's phase
  $group        this feature's group
  $db_id        this feature's internal database ID

=cut

# 'This is called when creating a feature from scratch.  It does not have
# an inherited coordinate system.

sub new {
  my $package = shift;
  my ($factory,
      $srcseq,
      $start,$stop,
      $method,$source,
      $score,$fstrand,$phase,
      $group,$db_id,$group_id,
      $tstart,$tstop) = @_;

  my $self = bless { },$package;
  ($start,$stop) = ($stop,$start) if defined($fstrand) and $fstrand eq '-';

  my $class =  $group ? $group->class : 'Sequence';

  @{$self}{qw(factory sourceseq start stop strand class)} =
    ($factory,$srcseq,$start,$stop,$fstrand,$class);

  # if the target start and stop are defined, then we use this information to create 
  # the reference sequence
  # THIS SHOULD BE BUILT INTO RELSEGMENT
  if (0 && $tstart ne '' && $tstop ne '') {
    if ($tstart < $tstop) {
      @{$self}{qw(ref refstart refstrand)} = ($group,$start - $tstart + 1,'+');
    } else {
      @{$self}{'start','stop'} = @{$self}{'stop','start'};
      @{$self}{qw(ref refstart refstrand)} = ($group,$tstop + $stop - 1,'-');
    }

  } else {
    @{$self}{qw(ref refstart refstrand)} = ($srcseq,1,'+');
  }

  @{$self}{qw(type fstrand score phase group db_id group_id absolute)} =
    (Bio::DB::GFF::Typename->new($method,$source),$fstrand,$score,$phase,
     $group,$db_id,$group_id,$factory->{absolute});

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
  if ($self->absolute) {
    return Bio::DB::GFF::RelSegment::_to_strand($self->{fstrand});
  }
  return $self->SUPER::strand || Bio::DB::GFF::RelSegment::_to_strand($self->{fstrand});
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

=head2 display_id

 Title   : display_id
 Usage   : $display_id = $f->display_id([$display_id])
 Function: get or set the feature display id
 Returns : a Bio::DB::GFF::Featname object
 Args    : a new display_id (optional)
 Status  : Public

This method is an alias for group().  It is provided for
Bio::SeqFeatureI compatibility.

=cut

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

*info         = \&group;
*display_id   = \&group;
*display_name = \&group;

=head2 target

 Title   : target
 Usage   : $target = $f->target([$new_target])
 Function: get or set the feature target
 Returns : a Bio::DB::GFF::Homol object
 Args    : a new group (optional)
 Status  : Public

This method works like group(), but only returns the group if it
implements the start() method.  This is typical for
similarity/assembly features, where the target encodes the start and
stop location of the alignment.

The returned object is of type Bio::DB::GFF::Homol, which is a
subclass of Bio::DB::GFF::Segment.

=cut


sub target {
    my $self = shift;
    my $group = $self->group or return;
    return unless $group->can('start');
    $group;
}

=head2 flatten_target

 Title   : flatten_target
 Usage   : $target = $f->flatten_target($f->target)
 Function: flatten a target object
 Returns : a string (GFF2), an array [GFF2.5] or an array ref [GFF3]
 Args    : a target object (required), GFF version (optional) 
 Status  : Public

This method flattens a target object into text for
GFF dumping.  If a second argument is provided, version-specific
vocabulary is used for the flattened target.

=cut

sub flatten_target {
    my $self = shift;
    my $t    = shift || return;
    my $v    = shift;

    return 0 unless $t->can('start');
    my $class = $t->class;
    my $name  = $t->name;
    my $start = $t->start;
    my $stop  = $t->stop;

    $v ||=2;
    if ( $v == 2.5 ) {
	
	print STDERR qq(Target "$class:$name"), "tstart $start", "tstop $stop\n";
	return (qq(Target "$class:$name"), "tstart $start", "tstop $stop");
    }
    elsif ( $v == 3 ) {
	return [Target=>"$name $start $stop"];
    }
    else {
	return qq(Target "$class:$name" $start $stop);
    }
}

# override parent a smidgeon so that setting the ref for top-level feature
# sets ref for all subfeatures
sub refseq {
  my $self   = shift;
  my $result = $self->SUPER::refseq(@_);
  if (@_) {
    my $newref = $self->SUPER::refseq;
    for my $sub ($self->get_SeqFeatures) {
      $sub->refseq(@_);
    }
  }
  $result;
}


=head2 hit

 Title   : hit
 Usage   : $hit = $f->hit([$new_hit])
 Function: get or set the feature hit
 Returns : a Bio::DB::GFF::Featname object
 Args    : a new group (optional)
 Status  : Public

This is the same as target(), for compatibility with
Bio::SeqFeature::SimilarityPair.

=cut

*hit = \&target;

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

sub id          { shift->{db_id}   }
sub primary_id  { shift->{db_id}   }

=head2 group_id

 Title   : group_id
 Usage   : $id = $f->group_id
 Function: get the feature ID
 Returns : a database identifier
 Args    : none
 Status  : Public

This method retrieves the database group identifier for the feature.
It cannot be changed.  Often the group identifier is more useful than
the feature identifier, since it is used to refer to a complex object
containing subparts.

=cut

sub group_id  { shift->{group_id}   }

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

  if (my $merged = $self->{merged_segs}) {
    $clone->{merged_segs} = { %$merged };
  }

  $clone;
}

=head2 compound

 Title   : compound
 Usage   : $flag = $f->compound([$newflag])
 Function: get or set the compound flag
 Returns : a boolean
 Args    : a new flag (optional)
 Status  : Public

This method gets or sets a flag indicated that the feature is not a
primary one from the database, but the result of aggregation.

=cut

sub compound  {
  my $self = shift;
  my $d    = $self->{compound};
  $self->{compound} = shift if @_;
  $d;
}

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feat = $feature->sub_SeqFeature([$method])
 Function: get subfeatures
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : a feature method (optional)
 Status  : Public

This method returns a list of any subfeatures that belong to the main
feature.  For those features that contain heterogeneous subfeatures,
you can retrieve a subset of the subfeatures by providing a method
name to filter on.

This method may also be called as segments() or get_SeqFeatures().

=cut

sub sub_SeqFeature {
  my $self = shift;
  my $type = shift;
  my $subfeat = $self->{subfeatures} or return;
  $self->sort_features;
  my @a;
  if ($type) {
    my $features = $subfeat->{lc $type} or return;
    @a = @{$features};
  } else {
    @a = map {@{$_}} values %{$subfeat};
  }
  return @a;
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

sub attach_seq { }


=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location 
	   of feature on sequence or parent feature  
 Returns : Bio::LocationI object
 Args    : none

=cut

sub location {
   my $self = shift;
   require Bio::Location::Split unless Bio::Location::Split->can('new');
   require Bio::Location::Simple unless Bio::Location::Simple->can('new');

   my $location;
   if (my @segments = $self->segments) {
       $location = Bio::Location::Split->new(-seq_id => $self->seq_id);
       foreach (@segments) {
          $location->add_sub_Location($_->location);
       }
   } else {
       $location = Bio::Location::Simple->new(-start  => $self->start,
					      -end    => $self->stop,
					      -strand => $self->strand,
					      -seq_id => $self->seq_id);
   }
   $location;
}

=head2 entire_seq

 Title   : entire_seq
 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns : a Bio::PrimarySeqI compliant object, or undef if there is no
           sequence attached
 Args    : none


=cut

sub entire_seq {
    my $self = shift;
    $self->factory->segment($self->sourceseq);
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

#'

sub merged_segments {
  my $self = shift;
  my $type = shift;
  $type ||= '';    # prevent uninitialized variable warnings

  my $truename = overload::StrVal($self);

  return @{$self->{merged_segs}{$type}} if exists $self->{merged_segs}{$type};
  my @segs = map  { $_->[0] } 
             sort { $a->[1] <=> $b->[1] ||
		    $a->[2] cmp $b->[2] }
             map  { [$_, $_->start, $_->type] } $self->sub_SeqFeature($type);

  # attempt to merge overlapping segments
  my @merged = ();
  for my $s (@segs) {
    my $previous = $merged[-1] if @merged;
    my ($pscore,$score) = (eval{$previous->score}||0,eval{$s->score}||0);
    if (defined($previous)
	&& $previous->stop+1 >= $s->start
	&& $pscore == $score
	&& $previous->method eq $s->method
       ) {
      if ($self->absolute && $self->strand < 0) {
	$previous->{start} = $s->{start};
      } else {
	$previous->{stop} = $s->{stop};
      }
      # fix up the target too
      my $g = $previous->{group};
      if ( ref($g) &&  $g->isa('Bio::DB::GFF::Homol')) {
	my $cg = $s->{group};
	$g->{stop} = $cg->{stop};
      }
    }
     elsif (defined($previous)
	    && $previous->start == $s->start
	    && $previous->stop  == $s->stop
	    && $previous->method eq $s->method
	   ) {
       next;
     }

  else {
      my $copy = $s->clone;
      push @merged,$copy;
    }
  }
  $self->{merged_segs}{$type} = \@merged;
  @merged;
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

=head2 attributes

 Title   : attributes
 Usage   : @attributes = $feature->attributes($name)
 Function: get the "attributes" on a particular feature
 Returns : an array of string
 Args    : feature ID
 Status  : public

Some GFF version 2 files use the groups column to store a series of
attribute/value pairs.  In this interpretation of GFF, the first such
pair is treated as the primary group for the feature; subsequent pairs
are treated as attributes.  Two attributes have special meaning:
"Note" is for backward compatibility and is used for unstructured text
remarks.  "Alias" is considered as a synonym for the feature name.

 @gene_names = $feature->attributes('Gene');
 @aliases    = $feature->attributes('Alias');

If no name is provided, then attributes() returns a flattened hash, of
attribute=E<gt>value pairs.  This lets you do:

  %attributes = $db->attributes;

=cut

sub attributes {
  my $self = shift;
  my $factory = $self->factory;
  defined(my $id = $self->id) or return;
  $factory->attributes($id,@_)
}


=head2 notes

 Title   : notes
 Usage   : @notes = $feature->notes
 Function: get the "notes" on a particular feature
 Returns : an array of string
 Args    : feature ID
 Status  : public

Some GFF version 2 files use the groups column to store various notes
and remarks.  Adaptors can elect to store the notes in the database,
or just ignore them.  For those adaptors that store the notes, the
notes() method will return them as a list.

=cut

sub notes {
  my $self = shift;
  $self->attributes('Note');
}

=head2 aliases

 Title   : aliases
 Usage   : @aliases = $feature->aliases
 Function: get the "aliases" on a particular feature
 Returns : an array of string
 Args    : feature ID
 Status  : public

This method will return a list of attributes of type 'Alias'.

=cut

sub aliases {
  my $self = shift;
  $self->attributes('Alias');
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

=head2 SeqFeatureI methods

The following Bio::SeqFeatureI methods are implemented:

primary_tag(), source_tag(), all_tags(), has_tag(), each_tag_value() [renamed get_tag_values()].

=cut

*primary_tag = \&method;
*source_tag  = \&source;
sub all_tags {
  my $self = shift;
  my %atts = $self->attributes;
  my @tags = keys %atts;

  # autogenerated methods
  #if (my $subfeat = $self->{subfeatures}) {
  #  push @tags,keys %$subfeat;
  #}

  @tags;
}
*get_all_tags = \&all_tags;

sub has_tag {
  my $self = shift;
  my $tag  = shift;
  my %att  = $self->attributes;
  my %tags = map {$_=>1} ( $self->all_tags );
  
  return $tags{$tag};
}

*each_tag_value = \&get_tag_values;

sub get_tag_values {
  my $self = shift;
  my $tag  = shift;
  return $self->$tag() if $CONSTANT_TAGS{$tag};
  
  my $atts = $self->attributes;
  return @{$atts->{$tag}} if $atts && $atts->{$tag};

  $tag = ucfirst $tag;
  return $self->$tag();  # try autogenerated tag
}

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  my $sub = $AUTOLOAD;
  my $self = $_[0];

  # ignore DESTROY calls
  return if $func_name eq 'DESTROY';

  # fetch subfeatures if func_name has an initial cap
#  return sort {$a->start <=> $b->start} $self->sub_SeqFeature($func_name) if $func_name =~ /^[A-Z]/;
  return $self->sub_SeqFeature($func_name) if $func_name =~ /^[A-Z]/;

  # error message of last resort
  $self->throw(qq(Can't locate object method "$func_name" via package "$pack"));
}#'

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
  my $shrink = shift;
  my $g = $self->{group};

  my $first = 0;
  my $tfirst = 0;
  if (my $subfeat = $self->{subfeatures}) {
    for my $list (values %$subfeat) {
      for my $feat (@$list) {
	# fix up our bounds to hold largest subfeature
	my($start,$stop,$strand) = $feat->adjust_bounds($shrink);

	if (defined($self->{fstrand})) {
	  $self->debug("Subfeature's strand ($strand) doesn't match parent strand ($self->{fstrand})\n") if $self->{fstrand} ne $strand;
	} else {
	  $self->{fstrand} = $strand;
	}

	my ($low,$high)  = $start < $stop ? ($start,$stop) : ($stop,$start);
	if ($shrink && !$first++) {
	  # first subfeature resets start & stop:
	  $self->{start} = $self->{fstrand} ne '-' ? $low : $high;
	  $self->{stop}  = $self->{fstrand} ne '-' ? $high : $low;
	} else {
	  if ($self->{fstrand} ne '-') {
	    $self->{start} = $low
	      if (!defined($self->{start})) || $low < $self->{start};
	    $self->{stop}  = $high
	      if (!defined($self->{stop}))  || $high  > $self->{stop};
	  } else {
	    $self->{start} = $high
	      if (!defined($self->{start})) || $high > $self->{start};
	    $self->{stop}  = $low
	      if (!defined($self->{stop}))  || $low  < $self->{stop};
	  }
	}

	# fix up endpoints of targets too (for homologies only)
	my $h = $feat->group;
	next unless $h && $h->isa('Bio::DB::GFF::Homol');
	next unless $g && $g->isa('Bio::DB::GFF::Homol');

	($start,$stop) = ($h->{start},$h->{stop});
	if ($shrink && !$tfirst++) {
	    $g->{start} = $start;
	    $g->{stop}  = $stop;
	} else {
	  if ($start <= $stop) {
	    $g->{start} = $start if (!defined($g->{start})) || $start < $g->{start};
	    $g->{stop}  = $stop  if (!defined($g->{stop}))  || $stop  > $g->{stop};
	  } else {
	    $g->{start} = $start if (!defined($g->{start})) || $start > $g->{start};
	    $g->{stop}  = $stop  if (!defined($g->{stop}))  || $stop  < $g->{stop};
	  }
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
      $subfeat->{$type} = [map { $_->[0] }
			   sort {$a->[1] <=> $b->[1] }
			   map { [$_,$_->start] }
			   @{$subfeat->{$type}}] if $strand > 0;
      $subfeat->{$type} = [map { $_->[0] }
			   sort {$b->[1] <=> $a->[1]}
			   map { [$_,$_->start] }
			   @{$subfeat->{$type}}] if $strand < 0;
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

sub name {
  my $self =shift;
  return $self->group || $self->SUPER::name;
}

=head2 gff_string

 Title   : gff_string
 Usage   : $string = $feature->gff_string
 Function: return GFF2 of GFF2.5 representation of feature
 Returns : a string
 Args    : none
 Status  : Public

=cut

sub gff_string {
  my $self = shift;
  my $version = $self->version;

  # gff3_string and gff_string are synonymous if the version is set to 3
  return $self->gff3_string(@_) if $version == 3;

  my ($start,$stop) = ($self->start,$self->stop);

  # the defined() tests prevent uninitialized variable warnings, when dealing with clone objects
  # whose endpoints may be undefined
  ($start,$stop) = ($stop,$start) if defined($start) && defined($stop) && $start > $stop;

  my ($class,$name) = ('','');
  my $strand = ('-','.','+')[$self->strand+1];

  my @group;

  if (my $t = $self->target) {
    push @group, $version == 2.5 ? $self->flatten_target($t,2.5) 
                                 : $self->flatten_target($t);
  }
  elsif (my $g = $self->group) {
    $class = $g->class || '';
    $name  = $g->name  || '';
    ($name =~ /\S\s\S/)?(push @group, "$class '$name'"):(push @group,"$class $name");
  }

  # add exhaustive list of attributes
  my $att = $self->attributes;
  for ( keys %$att ) {
      for my $v ( @{$att->{$_}} ) {     
	  $v = qq("$v") if $v=~ /\S\s+\S/;
	  push @group, qq($_ $v);
      }
  }

  my $group_field = join ' ; ',@group;
  my $ref = $self->refseq;
  my $n   = ref($ref) ? $ref->name : $ref;
  my $phase = $self->phase;
  $phase = '.' unless defined $phase;
  return join("\t",
	      $n,
	      $self->source,$self->method,
	      (defined $start ? $start : '.'),
	      (defined $stop  ? $stop  : '.'),
	      (defined $self->score ? $self->score : '.'),
	      (defined $strand ? $strand : '.'),
	      $phase,
	      $group_field);
}

=head2 gff3_string

 Title   : gff3_string
 Usage   : $string = $feature->gff3_string([$recurse])
 Function: return GFF3 representation of feature
 Returns : a string
 Args    : An optional flag, which if true, will cause the feature to recurse over
           subfeatures.
 Status  : Public

=cut

sub gff3_string {
  my $self = shift;
  my ($recurse,$parent) = @_;
  my ($start,$stop) = ($self->start,$self->stop);

  # the defined() tests prevent uninitialized variable warnings, when dealing with clone objects
  # whose endpoints may be undefined
  ($start,$stop) = ($stop,$start) if defined($start) && defined($stop) && $start > $stop;

  my $strand = ('-','.','+')[$self->strand+1];
  my $ref = $self->refseq;
  my $n   = ref($ref) ? $ref->name : $ref;
  my $phase = $self->phase;
  $phase = '.' unless defined $phase;

  my ($class,$name) = ('','');
  my @group;
  if (my $g = $self->group) {
    $class = $g->class || '';
    $name  = $g->name  || '';
    $name  = "$class:$name" if defined $class;
    push @group,[ID =>  $name] if !defined($parent) || $name ne $parent;
  }

  push @group,[Parent => $parent] if defined $parent && $parent ne '';

  if (my $t = $self->target) {
    $strand = '-' if $t->stop < $t->start;
    push @group, $self->flatten_target($t,3);
  }

  my @attributes = $self->attributes;
  while (@attributes) {
    push @group,[shift(@attributes),shift(@attributes)]
  }
  my $group_field = join ';',map {join '=',_escape($_->[0]),_escape($_->[1])} @group;
  my $string = join("\t",$n,$self->source,$self->method,$start||'.',$stop||'.',
                         $self->score||'.',$strand||'.',$phase,$group_field);
  $string .= "\n";
  if ($recurse) {
    foreach ($self->sub_SeqFeature) {
      $string .= $_->gff3_string(1,$name);
    }
  }
  $string;
}

=head2 version

 Title   : version
 Usage   : $feature->version()
 Function: get/set the GFF version to be returned by gff_string
 Returns : the GFF version (default is 2)
 Args    : the GFF version (2, 2.5 of 3)
 Status  : Public

=cut

sub version {
  my ($self, $version) = @_;
  $self->{version} = $version if $version;
  return $self->{version} || 2;
}


sub _escape {
  my $toencode = shift;
  $toencode    =~ s/([^a-zA-Z0-9_. :?^*\(\)\[\]@!-])/uc sprintf("%%%02x",ord($1))/eg;
  $toencode    =~ tr/ /+/;
  $toencode;
}

=head2 cmap_link()

 Title   : cmap_link
 Usage   : $link = $feature->cmap_link
 Function: returns a URL link to the corresponding feature in cmap
 Returns : a string
 Args    : none
 Status  : Public

If integrated cmap/gbrowse installation, it returns a link to the map otherwise
it returns a link to a feature search on the feature name.  See the cmap
documentation for more information.

This function is intended primarily to be used in gbrowse conf files. 
For example:

  link       = sub {my $self = shift; return $self->cmap_viewer_link(data_source);}

=cut


sub cmap_viewer_link {
  # Use ONLY if CMap is installed 
  my $self        = shift;
  my $data_source = shift;
  my $group_id    = $self->group_id;
  my $factory     = $self->factory; # aka adaptor

  my $link_str; 

  if ($factory->can("create_cmap_viewer_link")){
    $link_str = $factory->create_cmap_viewer_link(
        data_source => $data_source,
        group_id    => $group_id,
    );
  }
  my $name = $self->name();
  $link_str = '/cgi-bin/cmap/feature_search?features='
    . $name
    . '&search_field=feature_name&order_by=&data_source='
    . $data_source
    . '&submit=Submit'
    unless $link_str;

  return $link_str; 

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

