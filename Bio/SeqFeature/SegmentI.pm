package Bio::SeqFeature::SegmentI;

# $Id$
# A Bio::SeqFeature::CollectionI that is also a Bio::RelRangeI.

=head1 NAME

Bio::SeqFeature::SegmentI -- A segment of sequence that contains features.

=head1 SYNOPSIS

# get a Bio::SeqFeature::SegmentI somehow
# perhaps a Bio::SeqFeature::SimpleSegment

 use Bio::SeqFeature::SimpleSegment;
 my $segment =
   new Bio::SeqFeature::SimpleSegment(
     '-start' => 10,
     '-end' => 3000,
     '-strand' => 1,
     '-seq_id' => $seq_id
   );
 $segment->add_features( @featurelist );

 # Now the -range passed to features(..) will be relative to the
 # location [10,3000] on the forward strand.
 $segment->features(
    -range => new Bio::Location::Simple( '-start' => 1, '-end' => 300 ),
    -rangetype =>'overlaps'
 );

=head1 DESCRIPTION

The Bio::SeqFeature::SegmentI interface represents a segment of a
sequence that can contain features.  As such it implements the
L<Bio::RelRangeI> interface and the L<Bio::SeqFeature::CollectionI>
interface.  It uses RelRangeI vs RangeI to support switching between
relative and absolute positioning.

The only new method is segments(), which is an alias for features().

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.
Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...
use strict;
use vars qw( @ISA );

use Bio::RelRangeI;
use Bio::SeqFeature::CollectionI;
use Bio::DB::SegmentProviderI;
use Bio::LocallyIdentifiableI;
@ISA = qw( Bio::RelRangeI
           Bio::SeqFeature::CollectionI
           Bio::DB::SegmentProviderI
           Bio::LocallyIdentifiableI );

use vars '$VERSION';
$VERSION = '1.00';

use constant DEBUG_ADJUST_BOUNDS => 0;

=head2 unique_id (from Bio::LocallyIdentifiableI)

 Title   : unique_id
 Usage   : $id = $segment->unique_id( [$new_id] )
 Function: Getter/setter for the unique id for this segment
 Returns : the current (or former, if used as a set method) unique identifier
           (or undef)
 Args    : (optional) A new unique identifier, or "undef" if there is none
 Status  : Public

  This method will return undef if a unique identifier has not been
  set for this segment.  If the argument is the string "undef" then
  the unique_id will become undefined.  Note that the unique_id may
  not be changed (if, for instance, the implementing class does not
  allow unique_id changes).

=cut

sub unique_id {
  shift->throw_not_implemented( @_ );
}

#                   --Coders beware!--
# This pod is a modification of the pod for features() in
#   Bio/SeqFeature/CollectionI.pm
# , so changes must be kept in sync.
# Changes to this features() pod need to be copied to the following
# places, perhaps respecting existing modifications:
#   Bio/SeqFeature/SimpleSegment.pm [no  modifications]
#   Bio/SeqFeatureI.pm              [no  modifications]

=head2 features

 Title   : features
 Usage   : @features = $segment->features( %args );
           OR
           @features = $segment->features( @types );
 Returns : a list of L<Bio::SeqFeature::SegmentI> objects,
           OR
           (when the -iterator option is true) an L<Bio::SeqFeature::IteratorI>
           OR
           (when the -callback argument is given) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public

NOTE: This method is (almost) identical to the features() method
from L<Bio::SeqFeature::CollectionI> that it overrides.  The entire
documentation follows, but first a brief summary of the changes:
  * This method respects the absolute() flag of this SegmentI.  When
    absolute() is true the -baserange defaults to the abs_seq_id().
    When absolute() is false the -baserange defaults to $self, since
    this SegmentI isa L<Bio::RelRangeI>.
  * This method accepts the additional -absolute flag argument, which
    has the same effect as temporarily setting this SegmentI's
    absolute() flag (but is friendlier in a concurrent environment).
  * Features that are returned in relative mode (relative either to
    this SegmentI or to a given RangeI) will be returned with
    coordinates that are relative (and with their seq_id set to
    whatever it is that they are relative to, ie. whatever baserange
    ends up being).  Features that are returned in absolute mode will
    be returned with absolute coordinates.  The mode is determined by
    the -baserange and -absolute arguments and by the absolute() flag,
    in that precedence order.
  * If -rangetype is given but no -range is given then a special and
    strange thing happens: the method call is delegated to the
    parent_segment_provider.  If it is a SegmentI then its features()
    method will be called with all the same arguments but with *this*
    segment as the -range argument.  If the parent_segment_provider is
    a L<Bio::DB::SegmentProviderI> (but not a SegmentI) then the same
    thing will happen, but to the SegmentI returned by its
    get_collection() method with no arguments.  If the
    parent_segment_provider is null then no features will be returned.

This routine will retrieve features associated with this segment
object.  It can be used to return all features, or a subset based on
their type, location, or attributes.  Features that are returned in
relative mode (relative either to this SegmentI or to a given RangeI)
will be returned with coordinates that are relative.  Features that
are returned in absolute mode will be returned with absolute
coordinates.  The mode is determined by the -baserange and -absolute
arguments and by the absolute() flag, in that precedence order.

If ranges are specified using the -ranges argument, then these ranges
will be used to narrow the results, according to the specified
-rangetype and -strandtype arguments.

If no ranges are specified but the -rangetype argument is given then a
special and strange thing happens: the method call is delegated to the
parent_segment_provider.  If it is a SegmentI then its features()
method will be called with all the same arguments but with *this*
segment as the -range argument.  If the parent_segment_provider is a
L<Bio::DB::SegmentProviderI> (but not a SegmentI) then the same thing
will happen, but to the SegmentI returned by its get_collection()
method with no arguments.  If the parent_segment_provider is null then
no features will be returned.

If a -baserange is specified then unqualified ranges given with the
-ranges argument will be interpreted as relative to that baserange,
and qualified ranges will be re-relativized to the baserange.  If no
-baserange is given then a default will be provided that will depend
on the value of the -absolute argument or the absolute() flag.  If
-absolute is given and true or if absolute() is true then the default
baserange is the value returned by the abs_seq_id() method; if
( -absolute || absolute() ) is false then the default is this SegmentI
object ($self).  You may force absolute range interpretations by
giving a -baserange that is not a L<Bio::RangeI> (such as the string
'absolute', though any string will do the trick), by providing a true
value to the -absolute argument, or by setting the absolute() flag to
true.

-rangetype is one of:
   "overlaps"      return all features that overlap the range (default)
   "contains"      return features completely contained within the range
   "contained_in"  return features that completely contain the range

-strandmatch is one of:
   "strong"        ranges must have the same strand
                   (default ONLY when -strand is specified and non-zero)
   "weak"          ranges must have the same strand or no strand
   "ignore"        ignore strand information
                   (default unless -strand is specified and non-zero)

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of feature types (as if they
were given as -types => \@_).  In the named parameter form, the
arguments are a series of -name=E<gt>value pairs.  Note that the table
below is not exhaustive; implementations must support these but may
support other arguments as well (and are responsible for documenting the
difference).

  Argument       Description
  --------       ------------

  -type          A type name or an object of type L<Bio::SeqFeature::TypeI>
  -types         An array reference to multiple type names or TypeI objects

  -unique_id     A (string) unique_id.  See also -namespace.
  -unique_ids    An array reference to multiple unique_id values.

  -name          A (string) display_name or unique_id.  See also -namespace.
  -names         An array reference to multiple display_name/unique_id values.

  -namespace     A (string) namespace qualifier to help resolve the name/id(s)
  -class         same as -namespace

  -attributes    A hashref containing a set of attributes to match.  See
                 below.

  -baserange     A L<Bio::RangeI> object defining the range to which
                 the -range argument is relative.  The default
                 baserange depends on the value of the absolute()
                 flag.  If absolute() is true then the default is the
                 value of the abs_seq_id() method.  If absolute() is
                 false then the default is $self.  Note that the
                 baserange affects the sort order.  See also
                 -absolute.

  -absolute      If -absolute is given and true then all behavior will be as
                 if this SegmentI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SegmentI's
                 absolute() flag was set to false, even if it isn't.
                 Note that -baserange can still be given and can force
                 relativeness, and that takes precedence over -absolute.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch, -rangetype, and -baserange.
  -ranges        An array reference to multiple ranges.

  -rangetype     One of "overlaps", "contains", or "contained_in".  If no
                 range is given then a strange thing happens (it is
                 described above).

  -strandmatch   One of "strong", "weak", or "ignore".  Note that the
                 strand attribute of a given -range must be non-zero
                 for this to work (a 0/undef strand forces a 'weak'
                 strandmatch to become 'ignore' and cripples the
                 'strong' strandmatch).

  -iterator      Return a L<Bio::SeqFeature::IteratorI>

  -callback      A callback to invoke on each feature

  -sort          Return the features in order (of their start positions).
                 Note that if sorted() is true, then this argument is
                 redundant (the features will be returned in order
                 regardless).  If the baserange (see -baserange) has a
                 negative strand then the sort order will be reversed.

All plural arguments are interchangeable with their singular counterparts.

The -attributes argument is a hashref containing one or more
attributes to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple string matching, and multiple attributes
are ANDed together.  More complex filtering can be performed using the
-callback option (see below).

The -unique_ids argument is a reference to a list of strings.  Every
returned feature must have its unique_id value in this list or, if a
feature has no defined unique_id, then its display_name value in the
list if the list is provided.  A -unique_id argument is treated as a
single-element list of unique_ids.

The -names argument is a reference to a list of strings.  Every
returned feature must have its display_name or its unique_id value in this
list if the list is provided.  A -name argument is treated as a
single-element list of names.

If a -namespace is provided then names and ids (both queries and
targets) will be prepended with "$namespace:" as a bonus.  So
if you do features( -names => [ 'foo', 'bar' ], -namespace => 'ns' )
then any feature with the display_name or unique_id 'foo', 'ns:foo',
'bar', or 'ns:bar' will be returned.

If -iterator is true, then the method returns an object of type
Bio::SeqFeature::IteratorI.  Each call to next_seq() on this
object returns a Bio::SeqFeature::SegmentI object from this collection.

If -callback is passed a code reference, the code reference will be
invoked on each feature returned.  The code will be passed two
arguments consisting of the current feature and this SegmentI
object, and must return a true value. If the code returns a false
value, feature retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

-callback and -sort are mutually exclusive options.  If -sort is
defined, then -callback is ignored.  If you want to do a sorted
callback, set the sorted() flag of this collection to true.

If -sort or sorted() is true then the features will be returned in
order of the features' start positions.  This order will be reversed
if the baserange has a negative strand (remember that the default
baserange depends upon the value of the absolute() flag, but this may
be overridden by the -baserange argument).

Note that no guarantees are made by the SegmentI interface about
the order of the features, except when the sorted() flag is true or
when the -sort option is given to the features method.  Therefore
the implementation may choose to reorder the underlying data structure
to better accomodate -sorted feature requests as a result of a
features() call.  When this happens the SegmentI's sorted() flag
should be set to true, so that the client can detect that the -sorted
argument to features() is now irrelevant.

NOTE: the following methods all build on top of features(), and do not
need to be explicitly implemented.

    features_in_range()
    overlapping_features()
    contained_features()
    contained_in()
    get_feature_stream()
    get_feature_by_name()
    get_feature_by_id()
    get_feature_by_attribute()

=cut

sub features {
  shift->throw_not_implemented( @_ );
}

=head2 segments

 Title   : segments
 Usage   : @features = $segment->segments( %args );
           OR
           @features = $segment->segments( @types );
 Returns : a list of L<Bio::SeqFeature::SegmentI> objects,
           OR
           (when the -iterator option is true) an L<Bio::SeqFeature::IteratorI>
           OR
           (when the -callback argument is given) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public

  This method is a (glob ref) alias for features().

=cut

sub segments {
  shift->features( @_ );
}

#                   --Coders beware!--
# This pod is a modification of the pod for get_collection() in
#   Bio/DB/SegmentProviderI.pm
# , so changes must be kept in sync.
# Changes to this get_collection() pod need to be copied to the following
# places, perhaps respecting existing modifications:
#   Bio/SeqFeature/SimpleSegment.pm [no  modifications]
#   Bio/SeqFeatureI.pm              [no  modifications]

=head2 get_collection

 Title   : get_collection
 Usage   : my $segment = $segment->get_collection( %args );
           OR
           my $segment = $segment->get_collection( @types );
 Returns : A L<Bio::SeqFeature::SegmentI> object
 Args    : see below
 Status  : Public

NOTE: This method is (almost) identical to the get_collection() method
from L<Bio::DB::SegmentProviderI> that it overrides.  The entire
documentation follows, but first a brief summary of the changes:
  * This method respects the absolute() flag of this SegmentI.  When
    absolute() is true the -baserange defaults to the abs_seq_id().
    When absolute() is false the -baserange defaults to $self, since
    this SegmentI isa L<Bio::RelRangeI>.
  * This method accepts the additional -absolute flag argument, which
    has the same effect as temporarily setting this SegmentI's
    absolute() flag (but is friendlier in a concurrent environment).
  * SegmentIs returned will have as their seq_id() the baserange used
    here (if the baserange is absolute by any means then their seq_id()
    will be this SegmentI's abs_seq_id()).
  * This method will never throw the exception "Those features do not
    share a common sequence" because all features in this SegmentI
    already share a common sequence.

This routine will retrieve a L<Bio::SeqFeature::SegmentI> object based
on feature type, location or attributes.  The SeqFeatureI objects in
the returned SegmentI may or may not be newly instantiated by this
request.  They will have as their range the range searched, if any, or
the smallest range that encloses the returned features.  They will
have as their seq_id() the -baserange used here (if the baserange is
absolute by any means then their seq_id() will be this SegmentI's
abs_seq_id()).

If you make a modification to a feature you must call
update_collection with a collection that contains that feature to
ensure that the data provider is in sync with your change.  You may
not, however, assume that modifications to the feature do not
auto-sync (they might!).

If a range is specified using the -range argument then this range will
 be used to narrow the results, according to the specified -rangetype
 and -strandtype arguments.

-rangetype is one of:
   "overlaps"      return all features that overlap the range (default)
   "contains"      return features completely contained within the range
   "contained_in"  return features that completely contain the range

-strandmatch is one of:
   "strong"        ranges must have the same strand
   "weak"          ranges must have the same strand or no strand (default)
   "ignore"        ignore strand information

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of feature types (as if they
were given as -types => \@_).  In the named parameter form, the
arguments are a series of -name=E<gt>value pairs.  Note that the table
below is not exhaustive; implementations must support these but may
support other arguments as well (and are responsible for documenting the
difference).

  Argument       Description
  --------       ------------

  -type          A type name or an object of type L<Bio::SeqFeature::TypeI>
  -types         An array reference to multiple type names or TypeI objects

  -unique_id     A (string) unique_id.  See also -namespace.
  -unique_ids    An array reference to multiple unique_id values.

  -name          A (string) display_name or unique_id.  See also -namespace.
  -names         An array reference to multiple display_name/unique_id values.

  -namespace     A (string) namespace qualifier to help resolve the name/id(s)
  -class         same as -namespace

  -attributes    A hashref containing a set of attributes to match.  See
                 below.

  -baserange     A L<Bio::RangeI> object defining the range to which
                 the -range argument is relative.  The default
                 baserange depends on the value of the absolute()
                 flag.  If absolute() is true then the default is the
                 value of the abs_seq_id() method.  If absolute() is
                 false then the default is $self.  Note that the
                 baserange affects the sort order.  See also
                 -absolute.

  -absolute      If -absolute is given and true then all behavior will be as
                 if this SegmentI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SegmentI's
                 absolute() flag was set to false, even if it isn't.
                 Note that -baserange can still be given and can force
                 relativeness, and that takes precedence over -absolute.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch, -rangetype, and -baserange.
  -ranges        An array reference to multiple ranges.

  -rangetype     One of "overlaps", "contains", or "contained_in".

  -strandmatch   One of "strong", "weak", or "ignore".  Note that the strand
                 attribute of a given -range must be non-zero for this to work
                 (a 0/undef strand forces a 'weak' strandmatch to become
                 'ignore' and cripples the 'strong' strandmatch).

All plural arguments are interchangeable with their singular counterparts.

The -attributes argument is a hashref containing one or more
attributes to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple string matching, and multiple attributes
are ANDed together.

The -unique_ids argument is a reference to a list of strings.  Every
returned feature must have its unique_id value in this list or, if a
feature has no defined unique_id, then its display_name value in the
list if the list is provided.  A -unique_id argument is treated as a
single-element list of unique_ids.

The -names argument is a reference to a list of strings.  Every
returned feature must have its display_name or its unique_id value in this
list if the list is provided.  A -name argument is treated as a
single-element list of names.

If a -namespace is provided then names and ids (both queries and
targets) will be prepended with "$namespace:" as a bonus.  So
if you do features( -names => [ 'foo', 'bar' ], -namespace => 'ns' )
then any feature with the display_name or unique_id 'foo', 'ns:foo',
'bar', or 'ns:bar' will be returned.

=cut

sub get_collection {
  shift->throw_not_implemented( @_ );
}

=head2 adjust_bounds

 Title   : adjust_bounds
 Usage   : $segment->adjust_bounds( [ $shrink ], [ @sub_features ] );
 Function: Adjust the bounds of this feature to include all sub-features.
 Returns : nothing.
 Args    : [optional] a boolean value indicating whether or not to allow
           the bounds to shrink.
           [optional] a list of sub_features to adjust for.
 Status  : Public

  This method adjusts the boundaries of the feature to enclose all its
  sub_features (or, if a list of features is provided, then to enclose
  those).  If the $shrink argument is provided and true then the new
  boundaries will be the minimal bounds to include all sub_features.
  Otherwise the new bounds will never be smaller than the original
  bounds.  Note that, in either case, if the strand value of this
  feature is 0 when adjust_bounds() is called then it will be set to
  1.

=cut

## TODO: Test this! -Paul E
sub adjust_bounds {
  my $self = shift;
  my $shrink = shift;

  my $absolute = $self->absolute();
  if( $absolute ) {
    $self->absolute( 0 );
  }

  my @sub_features;
  my @all_features;
  if( defined( $shrink ) &&
      ref( $shrink ) &&
      $shrink->isa( 'Bio::SeqFeature::SegmentI' ) ) {
    @sub_features = ( $shrink, @_ );
    undef( $shrink );
  } else {
    @sub_features = @_;
  }
  unless( @sub_features ) {
    @all_features = $self->features();
    @sub_features = @all_features;
    ## TODO: REMOVE
    print STDERR "$self->adjust_bounds(..): using all subfeatures.\n" if DEBUG_ADJUST_BOUNDS;
  } else {
    ## TODO: REMOVE
    print STDERR "$self->adjust_bounds(..): \@sub_features are ( ", join( ", ", @sub_features ), " )\n" if DEBUG_ADJUST_BOUNDS;
  }
  # Don't adjust if there's no subfeatures.
  unless( @sub_features ) {
    ## TODO: REMOVE
    print STDERR "No subfeatures to adjust around.\n" if DEBUG_ADJUST_BOUNDS;
    # Done.  No change.
    if( $absolute ) {
      $self->absolute( 1 );
    }
    return;
  }

  my ( $self_abs_low,
       $self_abs_high,
       $self_abs_strand );
  unless( $shrink ) {
    ( $self_abs_low,
      $self_abs_high ) =
        ( $self->abs_low( 'plus' ),
          $self->abs_high( 'plus' ) );
  }
  my $totally_ignorant_of_self_bounds;
  unless( $self_abs_low && $self_abs_high &&
          ( defined( $self->start() ) || defined( $self->end() ) ) ) {
    undef $self_abs_low;
    undef $self_abs_high;
    $totally_ignorant_of_self_bounds = 1;
  }
  $self_abs_strand = $self->abs_strand();
  unless( $self->strand() ) {
    $self->strand( 1 );
  }
  my ( $sub_feature_abs_low,
       $sub_feature_abs_high );
  foreach my $sub_feature ( @sub_features ) {
    ## TODO: REMOVE
    print STDERR "                   adjusting for $sub_feature.\n" if DEBUG_ADJUST_BOUNDS;

    # As an added bonus we will clean up any sloppy features that
    # don't know where they are.  We do this by telling them that they
    # are here.  In this case, though, we should treat their coords as
    # having been absolute already.
    unless( defined $sub_feature->seq_id() ) {
      ## TODO: REMOVE
      print STDERR "                   -> setting its seq_id to $self.\n" if DEBUG_ADJUST_BOUNDS;
      ( $sub_feature_abs_low,
        $sub_feature_abs_high ) =
          ( $sub_feature->low( 'plus' ),
            $sub_feature->high( 'plus' ) );
      if( !$totally_ignorant_of_self_bounds &&
          defined( $self->start() ) ) {
        $sub_feature->start( 'plus',
                             ( $sub_feature_abs_low + 1 - $self->start( 'plus' ) ) );
        $sub_feature->end( 'plus',
                           ( $sub_feature_abs_high + 1 - $self->start( 'plus' ) ) );
      }
      if( $self->strand() < 0 ) {
        ## Now the strand is positive, related to me.
        $sub_feature->strand( 1 );
      }
      $sub_feature->seq_id( $self );
      ## TODO: REMOVE
      print STDERR "                   ->   (so now it is $sub_feature).\n" if DEBUG_ADJUST_BOUNDS;
    }

    ## TODO: REMOVE
    print STDERR "                   -> adjusting its bounds, with shrink = $shrink.\n" if DEBUG_ADJUST_BOUNDS;

    # All sub_features must recursively adjust their bounds.
    $sub_feature->adjust_bounds( $shrink );

    ## TODO: REMOVE
    print STDERR "                   ->   now it is $sub_feature.\n" if DEBUG_ADJUST_BOUNDS;

    # Skip any sub_features that are not on the same sequence.
    if( defined( $self->abs_seq_id() ) &&
        ( $sub_feature->abs_seq_id() ne $self->abs_seq_id() ) ) {
      ## TODO: Remove warning?
      warn "'$sub_feature', a subfeature of '$self', is defined on the sequence '", $sub_feature->abs_seq_id(), "', but '$self' is on the sequence '", $self->abs_seq_id(), "'.";
      next;
    }

    if( $totally_ignorant_of_self_bounds ) {
      ( $sub_feature_abs_low,
        $sub_feature_abs_high ) =
          ( $sub_feature->low( 'plus' ),
            $sub_feature->high( 'plus' ) );
    } else {
      ( $sub_feature_abs_low,
        $sub_feature_abs_high ) =
          ( $sub_feature->abs_low( 'plus' ),
            $sub_feature->abs_high( 'plus' ) );
    }
    if( !$self_abs_low ||
        ( $sub_feature_abs_low &&
          ( $sub_feature_abs_low < $self_abs_low ) )
      ) {
      ## TODO: REMOVE
      print STDERR "                   -> taking on its low value ($sub_feature_abs_low).\n" if DEBUG_ADJUST_BOUNDS;
      $self_abs_low = $sub_feature_abs_low;
      $totally_ignorant_of_self_bounds = 0;
    }
    if( !$self_abs_high ||
        ( $sub_feature_abs_high &&
          ( $sub_feature_abs_high > $self_abs_high ) ) ) {
      ## TODO: REMOVE
      print STDERR "                   -> taking on its high value ($sub_feature_abs_high).\n" if DEBUG_ADJUST_BOUNDS;
      $self_abs_high = $sub_feature_abs_high;
      $totally_ignorant_of_self_bounds = 0;
    }
  } # End foreach sub_feature.

  my ( $delta_low, $delta_high );
  unless( !$shrink && $self->abs_low( 'plus' ) && $self->abs_high( 'plus' ) &&
          ( defined( $self->start( 'plus' ) ) ||
            defined( $self->end( 'plus' ) ) )
        ) {
    $totally_ignorant_of_self_bounds = 1;
  }
  ## TODO: We should be able to remove the readjustment code below
  ## when all features are observers of their seq_feature parents and
  ## respond to the change in the parent's range by auto-adjusting..

  if( $totally_ignorant_of_self_bounds ) {
    ( $delta_low, $delta_high ) =
      ( $self_abs_low, $self_abs_high );
  } else {
    ( $delta_low, $delta_high ) =
      ( ( $self_abs_low - $self->abs_low( 'plus' ) ),
        ( $self_abs_high - $self->abs_high( 'plus' ) ) );
  }

  if( ( $delta_low == 0 ) && ( $delta_high == 0 ) ) {
    ## TODO: REMOVE
    print STDERR "                   no change.\n" if DEBUG_ADJUST_BOUNDS;
    # Done.  No change.
    if( $absolute ) {
      $self->absolute( 1 );
    }
    return;
  }

  ## TODO: REMOVE
  print STDERR "delta_low is $delta_low, delta_high is $delta_high.\n" if DEBUG_ADJUST_BOUNDS;

  # We're going to have to go through all features that we contain,
  # plus those passed in as arguments (which may or may not be
  # subfeatures, really, despite the fact that we store it in an array
  # called @sub_features, just to be confusing).
  unless( @all_features ) {
    @all_features = $self->features();
    my @new_features;
    foreach my $sub_feature ( @sub_features ) {
      unless( grep { $sub_feature == $_ } @all_features ) {
        push( @new_features, $sub_feature );
      }
    }
    if( @new_features ) {
      push( @all_features, @new_features );
    }
  }

  # Now change each sub_feature's bounds by the deltas.
  my ( $sub_feature_low,
       $sub_feature_high,
       $sub_feature_strand,
       $sub_feature_ancestor,
       $last_sub_feature_ancestor );
  foreach my $sub_feature ( @all_features ) {
    ## TODO: REMOVE
    print STDERR "                   adjusting $sub_feature for change to self.\n" if DEBUG_ADJUST_BOUNDS;
    # Skip any sub_features that are not physically contained herein.
    $last_sub_feature_ancestor = $sub_feature;
    $sub_feature_ancestor = $sub_feature->seq_id();
    while( ref( $sub_feature_ancestor ) &&
           $sub_feature_ancestor->isa( 'Bio::SeqFeature::SegmentI' ) &&
           !( $sub_feature_ancestor == $self ) ) {
      $last_sub_feature_ancestor = $sub_feature_ancestor;
      $sub_feature_ancestor = $sub_feature_ancestor->seq_id();
    }
    ## TODO: Put back.
    #next unless( $sub_feature_ancestor == $self );
    ## TODO: REMOVE.
    unless( $sub_feature_ancestor == $self ) {
      warn "EEgad, man, subfeature $sub_feature doesn't have self as an ancestor.";
      next;
    }

    ( $sub_feature_low,
      $sub_feature_high,
      $sub_feature_strand ) =
        ( $sub_feature->low( 'plus' ),
          $sub_feature->high( 'plus' ),
          $sub_feature->strand() );
    next unless $delta_low;
    ## TODO: REMOVE
    print STDERR "                   readjusting $sub_feature by ( $delta_low )\n" if DEBUG_ADJUST_BOUNDS;
    if( $totally_ignorant_of_self_bounds ) {
      ## If we had undefined bounds originally (or were otherwise
      ## ignorant of our bounds) then delta_low is effectively
      ## offset from 0, not 1.  We compensate by adding 1 here.
      $sub_feature->start( 'plus', ( $sub_feature_low + 1 - $delta_low ) );
      $sub_feature->end( 'plus', ( $sub_feature_high + 1 - $delta_low ) );
    } else {
      $sub_feature->start( 'plus', ( $sub_feature_low - $delta_low ) );
      $sub_feature->end( 'plus', ( $sub_feature_high - $delta_low ) );
    }
    $sub_feature->ensure_orientation();
    # Just in case its ancestor ceases to be == us after we adjust
    # ourself, we'll force it to be...
    $last_sub_feature_ancestor->seq_id( $self );
    ## TODO: REMOVE
    print STDERR "$sub_feature bounds is now ", $sub_feature->Bio::RelRangeI::toString(), "\n" if DEBUG_ADJUST_BOUNDS;
  } # End foreach sub_feature, rebound it so that when we change our
    # bounds, their absolute bounds remain the same.

  # Now, finally, change our own bounds.
  if( $totally_ignorant_of_self_bounds ) {
    ## If we're ignorant of our bounds then the delta values are our
    ## new absolute positions.  We might need to relativize the values
    ## before setting them.
    my $seq_id = $self->seq_id();
    if( defined( $seq_id ) &&
        ref( $seq_id ) &&
        $seq_id->isa( 'Bio::RangeI' ) 
      ) {
      # Relativize.
      my $seq_id_low;
      if( $seq_id->isa( 'Bio::RelRangeI' ) ) {
        $seq_id_low = $seq_id->low( 'plus' );
      } else {
        $seq_id_low = $seq_id->start();
      }
      $self->start( 'plus', ( $seq_id_low + $delta_low - 1 ) );
      $self->end( 'plus', ( $seq_id_low + $delta_high - 1 ) );
    } else {
      $self->start( 'plus', $delta_low );
      $self->end( 'plus', $delta_high );
    }
  } else {
    $self->start( 'plus', $self->start() + $delta_low );
    $self->end( 'plus', $self->end() + $delta_high );
  }
  $self->ensure_orientation();
  ## TODO: REMOVE
  print STDERR "$self bounds is now ", $self->Bio::RelRangeI::toString(), "\n" if DEBUG_ADJUST_BOUNDS;

  if( $absolute ) {
    $self->absolute( 1 );
  }
  return;
} # adjust_bounds(..)

## TODO: I dunno what.  Bio::DB::GFF::Browser expects this field to exist.
sub class {
  return;
}

1;

__END__
