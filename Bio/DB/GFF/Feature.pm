# 

=head1 NAME

Bio::DB::GFF::Feature -- A relative segment identified by a feature type

=head1 SYNOPSIS

See L<Bio::DB::GFF>.

=head1 DESCRIPTION

Bio::DB::GFF::Feature is a stretch of sequence that corresponding to a
single annotation in a GFF database.  It inherits from
Bio::DB::GFF::Segment, and so has all the support for relative
addressing of this class and its ancestors.  It also inherits from
Bio::SeqFeatureI and so has the familiar start(), stop(),
primary_tag() and location() methods (it implements Bio::LocationI
too, if needed).

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
Bio::DB::GFF::Segment-E<gt>features() method.

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

use Bio::SeqFeature::Generic;
use Bio::DB::GFF::Segment;

use vars qw( $VERSION @ISA $AUTOLOAD );
@ISA = qw( Bio::SeqFeature::Generic
           Bio::DB::GFF::Segment );

$VERSION = '0.61';

use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Featname;
use Bio::DB::GFF::Typename;
use Bio::DB::GFF::Homol;
use Bio::LocationI;

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
Segment object, and should inherit that object's coordinate system.

The 14 arguments are positional (sorry):

  $parent       a Bio::DB::GFF::Segment object (or descendent)
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
  $group_id     this feature's internal group database ID
  $tstart       this feature's target start
  $tstop        this feature's target stop

tstart and tstop aren't used for anything at the moment, since the
information is embedded in the group object.

=cut

# this is called for a feature that is attached to a parent sequence,
# in which case it inherits its coordinate reference system and strandedness
sub new_from_parent {
  my $pack   = shift;
  my ($parent,
      $start,$stop,
      $method,$source,$score,
      $fstrand,$phase,
      $group,$db_id,$group_id,
      $tstart,$tstop) = @_;

  if( $start > $stop ) {
    $fstrand = '-';
  }
  ($start,$stop) = ($stop,$start) if defined($fstrand) and $fstrand eq '-';
  unless( defined( $fstrand ) ) {
    $fstrand = 0;
  }
  my $class = $group ? $group->class : $parent->class;

  # The given parent is a Segment, but we want the factory aka
  # parent_segment_provider to be the GFF object at the top (grandad,
  # probably).
  my $seq_id = $parent;
  while( $parent && ref( $parent ) && !$parent->isa( 'Bio::DB::GFF' ) ) {
    last unless( $parent->can( 'factory' ) );
    $parent = $parent->factory();
  }

  my %args =
    ( '-seq_id' => $seq_id,
      '-start' => $seq_id->abs2rel( $start ),
      '-end' => $seq_id->abs2rel( $stop ),
      '-strand' => $seq_id->abs2rel_strand( $fstrand ),
      '-absolute' => $parent->absolute(),
      '-type' => Bio::DB::GFF::Typename->new( $method, $source ),
      '-parent' => $parent,
      '-orientation_policy' => 'dependent'
    );
  my $self = $pack->SUPER::new( %args );
  $self->score( $score ) if defined( $score );
  $self->phase( $phase ) if defined( $phase );
  $self->group( $group ) if defined( $group );
  $self->class( $class ) if defined( $class );
  $self->{ '_db_id' } = $db_id;
  $self->{ '_group_id' } = $group_id;
  $self->{ 'merged_seqs' } = undef;
  # Default to sorted
  $self->sorted( 1 );

  return $self;
} # new_from_parent(..)

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
without reference to a Segment object, and should therefore use its
default coordinate system (relative to itself).

The 14 arguments are positional:

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
  $group_id     this feature's internal group database ID
  $tstart       this feature's target start
  $tstop        this feature's target stop

tstart and tstop aren't used for anything at the moment, since the
information is embedded in the group object.

=cut

# 'This is called when creating a feature from scratch.  It does not have
# an inherited coordinate system.
sub new {
  my $pack = shift;
  my ($factory,
      $srcseq,
      $start,$stop,
      $method,$source,
      $score,$fstrand,$phase,
      $group,$db_id,$group_id,
      $tstart,$tstop) = @_;

  if( $start > $stop ) {
    $fstrand = '-';
  }
  ($start,$stop) = ($stop,$start) if defined($fstrand) and $fstrand eq '-';
  unless( defined( $fstrand ) ) {
    $fstrand = 0;
  }
  my $class =  $group ? $group->class : 'Sequence';

  my %args =
    ( '-seq_id' => $srcseq,
      '-start' => $start,
      '-end' => $stop,
      '-strand' => $fstrand,
      '-absolute' => ( defined( $factory ) ? $factory->absolute() : undef ),
      '-type' => Bio::DB::GFF::Typename->new( $method, $source ),
      '-parent' => $factory,
      '-orientation_policy' => 'dependent'
    );
  my $self = $pack->SUPER::new( %args );
  $self->{ '_db_id' } = $db_id;
  $self->{ '_group_id' } = $group_id;
  $self->{ 'merged_seqs' } = undef;
  $self->score( $score ) if defined( $score );
  $self->group( $group ) if defined( $group );
  $self->phase( $phase ) if defined( $phase );
  $self->class( $class ) if defined( $class );

  return $self;
} # new(..)

=head2 new_from_feature

 Title   : new_from_feature
 Usage   : my $new_feature =
             Bio::DB::GFF::Feature->new_from_feature( $copy_from );
 Function: Create a new Bio::DB::GFF::Feature object by copying
           values from another Feature object.
 Returns : A new L<Bio::DB::GFF::Feature> object
 Args    : Another L<Bio::DB::GFF::Feature> object
 Status  : Protected

  This is a special copy constructor.  It forces the new feature into
  the L<Bio::DB::GFF::Feature> package, regardless of the package that
  it is called from.  This causes subclass-specfic information to be
  dropped.

  This also does not copy into the new feature the features held in
  the existing feature.  If you would like the new feature to hold the
  same features you must explicitly add them, like so:
    $new_feature->add_features( $copy_from->features() );

  As a special bonus you may also pass an existing hash and it will be the
  blessed an anointed object that is returned, like so:
    $new_feature =
      Bio::DB::GFF::Feature->new_from_feature(
        $copy_from,
        $new_feature
      );

=cut

sub new_from_feature {
  my $pack = shift; # ignored
  my $feature = shift || $pack;
  my $new_feature = shift;
  $new_feature =
    Bio::SeqFeature::Generic->new_from_feature(
      $feature,
      $new_feature
    );
  $new_feature =
    Bio::DB::GFF::Segment->new_from_segment(
      $feature,
      $new_feature
    );
  @{ $new_feature }{ qw( _db_id _group_id merged_seqs ) } =
    @{ $feature }{ qw( _db_id _group_id merged_seqs ) };

  return bless $new_feature, __PACKAGE__;
} # new_from_feature(..)

=head2 new_from_segment

 Title   : new_from_segment
 Usage   : $s = Bio::DB::GFF::Feature->new_from_segment( $copy_from )
 Function: create a new L<Bio::DB::GFF::Segment>
 Returns : A new L<Bio::DB::GFF::Segment> object
 Args    : Another L<Bio::DB::GFF::Segment> object
 Status  : Protected

  This constructor is used internally by the subseq() method.  It forces
  the new segment into the L<Bio::DB::GFF::Segment> package, regardless
  of the package that it is called from.  This causes subclass-specific
  information, such as feature types, to be dropped when a subsequence
  is created.

  This also does not copy into the new segment the features held in
  the existing segment.  If you would like the new segment to hold the
  same features you must explicitly add them, like so:
    $new_segment->add_features( $copy_from->features() );

  As a special bonus you may also pass an existing hash and it will be the
  blessed an anointed object that is returned, like so:
    $new_segment =
      Bio::DB::GFF::Segment->new_from_segment(
        $copy_from,
        $new_segment
      );

  This delegates explicitly to the L<Bio::DB::GFF::Segment> superclass
  method of the same name.

=cut

## This method is here because a method of the same name exists in
## both superclasses, and we want to use the Bio::DB::GFF::Segment
## one..
sub new_from_segment {
  shift->Bio::DB::GFF::Segment::new_from_segment( @_ );
} # new_from_segment(..)

## These methods are here because a method of the same name exists in
## both superclasses, and we want to use the Bio::DB::GFF::Segment
## one..
sub seq {
  shift->Bio::DB::GFF::Segment::dna( @_ );
} # seq(..)
sub dna {
  shift->Bio::DB::GFF::Segment::dna( @_ );
} # dna(..)
sub protein {
  shift->Bio::DB::GFF::Segment::dna( @_ );
} # protein(..)

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
  my $type = $self->type();
  if( !defined( $type ) ) {
    if( @_ ) {
      $type = Bio::DB::GFF::Typename->new();
    } else {
      return;
    }
  } 
  if( $type->can( 'method' ) ) {
    my $d = $type->method();
    $type->method( shift ) if @_;
    $d;
  }
} # method

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
  my $type = $self->type;
  if( !defined( $type ) ) {
    if( @_ ) {
      $type = Bio::DB::GFF::Typename->new();
    } else {
      return;
    }
  } 
  if( $type->can( 'source' ) ) {
    my $d = $type->source();
    $type->source( shift ) if @_;
    $d;
  }
} # source

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
  shift->frame( @_ );
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

sub group {
  my $self = shift;
  my $d = $self->{ '_gsf_tag_hash' }->{'group'};
  $self->{ '_gsf_tag_hash' }->{'group'} = shift if @_;
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

*info         = \&group;
*display_id   = \&group;
*display_name = \&group;

=head2 target

 Title   : target
 Usage   : $target = $f->target([$new_target])
 Function: get or set the feature target
 Returns : a Bio::DB::GFF::Featname object
 Args    : a new group (optional)
 Status  : Public

This method works like group(), but only returns the group if it
implements the start() method.  This is typical for
similarity/assembly features, where the target encodes the start and stop
location of the alignment.

=cut

sub target {
  my $self = shift;
  my $group = $self->group or return;
  return unless $group->can('start');
  $group;
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

sub id        { shift->{ '_db_id' }   }

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

sub group_id  { shift->{ '_group_id' }   }

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
} # clone(..)

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

sub compound {
  my $self = shift;
  my $d = $self->{ '_gsf_tag_hash' }->{'compound'};
  $self->{ '_gsf_tag_hash' }->{'compound'} = shift if @_;
  $d;
}

=head2 add_subfeature

 Title   : add_subfeature
 Usage   : $feature->add_subfeature($feature)
 Function: add a subfeature to the feature
 Returns : nothing
 Args    : a Bio::DB::GFF::Feature object
 Status  : Public

An alias for add_features(..)

=cut

sub add_subfeature {
  shift->add_features( @_ );
} # add_subfeature(..)

### TODO: ERE I AM
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

  my $truename = overload::StrVal($self);

  return @{$self->{merged_segs}{$type}} if exists $self->{merged_segs}{$type};
  my @segs = map  { $_->[0] } 
             sort { $a->[1] <=> $b->[1] ||
		    $a->[2] cmp $b->[2] }
             map  { [$_, $_->low( 'plus' ), $_->type] } $self->sub_SeqFeature($type);

  # attempt to merge overlapping segments
  my @merged = ();
  for my $s (@segs) {
    my $previous = $merged[-1] if @merged;
    my ($pscore,$score) = (eval{$previous->score}||0,eval{$s->score}||0);
    if (defined($previous) 
	&& $previous->high( 'plus' )+1 >= $s->low( 'plus' )
	&& $pscore == $score
       ) {
      if ($self->absolute && $self->strand < 0) {
	$previous->low( 'plus', $s->low( 'plus' ) );
      } else {
	$previous->high( 'plus', $s->high( 'plus' ) );
      }
      # fix up the target too
      my $g = $previous->{group};
      if ( ref($g) &&  $g->isa('Bio::DB::GFF::Homol')) {
	my $cg = $s->{group};
	$g->end( 'plus', $cg->end( 'plus' ) );
      }
    } elsif (defined($previous) 
	     && $previous->low( 'plus' ) == $s->low( 'plus' ) 
	     && $previous->low( 'plus' ) == $s->low( 'plus' )) {
      next;
    } else {
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

An alias for types()

=cut

sub sub_types {
  return shift->types( @_ );
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

=head2 has_tag

 Title   : has_tag
 Usage   : $tag_exists = $self->has_tag('some_tag')
 Function: 
 Returns : TRUE if the specified tag exists, and FALSE otherwise
 Args    : a tag name
 Status  : Public

=cut

sub has_tag {
  my $self = shift;
  my ( $tag ) = @_;
  return $self->SUPER::has_tag( $tag ) || scalar( $self->attributes( $tag ) );
} # has_tag(..)

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : @values = $self->get_tag_values( 'some_tag' );
 Function: 
 Returns : An array comprising the values of the specified tag.
 Args    :


=cut

sub get_tag_values {
  my $self = shift;
  my ( $tag ) = @_;
  my @tag_values =
    ( $self->SUPER::has_tag( $tag ) ?
      $self->SUPER::get_tag_values( $tag ) :
      () );
  push( @tag_values, $self->attributes( $tag ) );
  return @tag_values;
} # get_tag_values(..)

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : @tags = $feat->get_all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none
 Status  : Public

=cut

sub get_all_tags {
  my $self = shift;
  my @tags = $self->SUPER::get_all_tags();
  my @all_attributes = $self->attributes();
  my %all_attributes_hash = @all_attributes;
  push( @tags, keys( %all_attributes_hash ) );
  return @tags;
} # get_all_tags()

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

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  my $sub = $AUTOLOAD;
  my $self = $_[0];

  # ignore DESTROY calls
  return if $func_name eq 'DESTROY';

  # fetch subfeatures if func_name has an initial cap
#  return sort {$a->low( 'plus' ) <=> $b->low( 'plus' )} $self->sub_SeqFeature($func_name) if $func_name =~ /^[A-Z]/;
  return $self->sub_SeqFeature($func_name) if $func_name =~ /^[A-Z]/;

  # error message of last resort
  $self->throw( "Can't locate object method \"$func_name\" via package \"$pack\"" );
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

This is an alias for sorted( 1 );

=cut

# sort features
sub sort_features {
  shift->sorted( 1 );
} # sort_features(..)

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
  ## TODO: REMOVE?
  #return $self->toRelRangeString( 'both', 'plus' );
#  my $type = $self->type;
#  my $name = $self->group;
#  return "$type($name)" if $name;
#  return $type;
  my $type = $self->method;
  my $id   = $self->group || 'unidentified';
  return join '/',$id,$type,$self->SUPER::asString;
}

sub name {
  my $self =shift;
  return $self->group() || $self->SUPER::display_name();
}

sub gff_string {
  my $self = shift;
  my ($start,$stop) = ($self->low( 'plus' ),$self->high( 'plus' ));

  # the defined() tests prevent uninitialized variable warnings, when dealing with clone objects
  # whose endpoints may be undefined
  ($start,$stop) = ($stop,$start) if defined($start) && defined($stop) && $start > $stop;

  my ($class,$name) = ('','');
  my @group;
  if (my $t = $self->target) {
    my $class = $t->class;
    my $name  = $t->name;
    my $start = $t->low( 'plus' );
    my $stop  = $t->high( 'plus' );
    push @group,qq(Target "$class:$name" $start $stop);
  }

  elsif (my $g = $self->group) {
    $class = $g->class || '';
    $name  = $g->name  || '';
    push @group,"$class $name";
  }
  push @group,map {qq(Note "$_")} $self->notes;

  my $group_field = join ' ; ',@group;
  my $strand = ('-','.','+')[$self->strand+1];
  my $ref = $self->refseq;
  my $n   = ref($ref) ? $ref->name : $ref;
  my $phase = $self->phase;
  $phase = '.' unless defined $phase;
  return join("\t",$n,$self->source,$self->method,$start||'.',$stop||'.',$self->score||'.',$strand||'.',$phase,$group_field);
}

=head1 A Note About Similarities

The current default aggregator for GFF "similarity" features creates a
composite Bio::DB::GFF::Feature object of type "gapped_alignment".
The target() method for the feature as a whole will return a
Segment object that is as long as the extremes of the similarity
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

L<bioperl>, L<Bio::DB::GFF>, L<Bio::DB::Segment>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

