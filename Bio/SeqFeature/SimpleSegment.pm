package Bio::SeqFeature::SimpleSegment;

# $Id $
# A simple implementation of a Bio::SeqFeature::CollectionI that is
# also a Bio::RelRangeI.  Implemented as a
# Bio::SeqFeature::SimpleCollection that is also a Bio::RelRange.

=head1 NAME

Bio::SeqFeature::SimpleSegment -- A segment of sequence that contains features.

=head1 SYNOPSIS

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

Bio::SeqFeature::SimpleSegment is a simple implementation of the
L<Bio::SeqFeature::SegmentI> interface, which represents a segment of
a sequence that can contain features.  This implementation extends
L<Bio::SeqFeature::SimpleCollection> and L<Bio::RelRange>, with some
extra code to support the enhancements to features() and
get_collection() required by the interface.

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

use Bio::RelRange;
use Bio::SeqFeature::SimpleCollection;
use Bio::SeqFeature::SegmentI;
@ISA = qw( Bio::RelRange
           Bio::SeqFeature::SimpleCollection
           Bio::SeqFeature::SegmentI );

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'

use vars '$VERSION';
$VERSION = '1.00';

=head2 new

  Title   : new
  Usage   : $range = Bio::SeqFeature::SimpleSegment->new(
                       -seq_id => $a_range_or_a_sequence_id,
                       -start => 100,
                       -end => 200,
                       -strand => +1,
                       -absolute => 1,
                       -features => $feature_list_ref,
                       -type => $a_type_or_string,
                       -parent => undef
                     );
  Function: generates a new Bio::SeqFeature::SimpleSegment object
  Returns : a new Bio::SeqFeature::SimpleSegment object
  Args    : two of (-start, -end, -length) - the third is calculated
          : -strand (defaults to 0)
          : -seq_id (a L<Bio::RangeI> or a (string) id of a sequence)
          : -absolute (defaults to 0)
          : -features (an optional array ref of L<Bio::SeqFeatureI> objects)
          : -type     (an optional L<Bio::SeqFeature::TypeI> or string)
          : -parent   (an optional L<Bio::DB::SegmentProviderI> parent object)

=cut

sub new {
  my ( $caller, @args ) = @_;

  # Rely on RelRange's new(..) to do most of the work.
  my $self = $caller->Bio::RelRange::new( @args );
  my ( $features, $type, $parent ) =
    $self->_rearrange( [ qw( FEATURES TYPE PARENT ) ], @args );
  if( $features ) {
    foreach my $feature ( @$features ) {
      next unless( ref( $feature ) && $feature->isa( "Bio::SeqFeatureI" ) );
      unless( $self->_insert_feature( $feature ) ) {
        $self->throw( "duplicate feature: $feature" );
      }
    }
  }
  $self->type( $type ) if defined( $type );
  $self->parent_segment_provider( $parent ) if defined( $parent );
  return $self;
} # new(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 features

 Title   : features
 Usage   : @features = $segment->features( %args );
           OR
           @features = $segment->features( @types );
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           (when the -iterator option is true) an L<Bio::SeqFeature::IteratorI>
           OR
           (when the -callback argument is given) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public

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

  -rangetype     One of "overlaps", "contains", or "contained_in".  If no
                 range is given then a strange thing happens (it is
                 described above).

  -strandmatch   One of "strong", "weak", or "ignore".  Note that the strand
                 attribute of a given -range must be non-zero for this to work
                 (a 0/undef strand forces a 'weak' strandmatch to become
                 'ignore' and cripples the 'strong' strandmatch).

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
object returns a Bio::SeqFeatureI object from this collection.

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
  my $self = shift;
  my ( $types, $absolute, $baserange, $ranges, $rangetype, $iterator, $callback );
  my %args;
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $types, $absolute, $baserange, $ranges, $rangetype, $iterator, $callback ) =
      rearrange(
        [ [ qw( TYPE TYPES ) ],
          [ qw( ABSOLUTE ABS ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ],
          [ qw( RANGETYPE RANGE_TYPE ) ],
          [ qw( ITERATOR STREAM ) ],
          [ qw( CALLBACK CALL_BACK ) ]
        ],
        @_
      );
    %args = %_;
  } else {
    ## Types.
    $types = \@_;
    %args = { -types => $types };
  }
  # If -baserange is given but is not a RangeI then the
  # -baserange argument to the superclass should be
  # $self->abs_seq_id().  If -baserange is not given then the argument
  # to the superclass method should be $self unless -absolute is given
  # and true or absolute() is true.
  if( not defined( $baserange ) ) {
    if( $absolute || $self->absolute() ) {
      $baserange = $self->abs_seq_id();
    } else {
      $baserange = $self;
    }
  } elsif( defined( $baserange ) &&
           ( not ref( $baserange ) || not $baserange->isa( 'Bio::RangeI' ) )
         ) {
    $baserange = $self->abs_seq_id();
  }
  $args{ '-baserange' } = $baserange;
  # Oh yeah and make sure it's the only baserange argument..
  $args{ '-base_range' } = $args{ '-baselocation' } =
    $args{ '-base_location' } = $args{ '-base' } =
      $args{ 'baserange' } = $args{ 'base_range' } =
        $args{ 'baselocation' } = $args{ 'base_location' } =
          $args{ 'base' } = undef;

  # If $rangetype is given but not $ranges then we delegate up to daddy.
  if( defined( $rangetype ) && not defined( $ranges ) ) {
    my $daddy = $self->parent_segment_provider();
    unless( defined $daddy ) {
      # If there's no parent then we must return the info that no
      # features match the request.
      if( $iterator ) {
        return Bio::SeqFeature::SimpleIterator->new();
      } elsif( $callback ) {
        return 1; # All okay.
      } else {
        return undef; # Empty list.
      }
    }
    my $daddy_segment;
    if( $daddy->isa( 'Bio::SeqFeature::SegmentI' ) ) {
      $daddy_segment = $daddy;
    } else {
      $daddy_segment = $daddy->get_collection();
    }
    $args{ '-range' } = $self;
    # PS we already know that it's the only range argument
    return $daddy_segment->( %args ); 
  } # End if $rangetype is defined but no $ranges are given.

  # We also have promised that the returned features will be relative to whatever
  # baserange ends up being.
  # There's 3 options for this: by default features() returns a
  # list. if -iterator is given then it returns an iterator.  if -callback
  # is given then some method is called.  We need to deal with each option
  # separately.
  if( $iterator ) {
    # Wrap the returned iterator in our own..
    my $super_iterator = $self->SUPER::features( %args );
    ## SimpleSegmentIteratorWrapper is defined in this file, below.
    return Bio::SeqFeature::SimpleSegmentIteratorWrapper->new(
      $baserange,
      $iterator
    );
  } elsif( $callback ) {
    # Wrap the callback in our own..
    my $new_callback = $self->_create_callback_wrapper( $baserange, $callback );
    $args{ '-callback' } = $new_callback;
    # Oh and make sure it's the only callback argument..
    $args{ '-call_back' } = $args{ 'callback' } = $args{ 'call_back' } = undef;
    return $self->SUPER::features( %args );
  } else {
    # post-process the list.
    return map
      { $self->_relativize_feature( $baserange, $_ ) }
        $self->SUPER::features( %args );
  } # End if $iterator .. or $callback .. or else.
} # features(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 get_collection

 Title   : get_collection
 Usage   : my $segment = $segment->get_collection( %args );
           OR
           my $segment = $segment->get_collection( @types );
 Returns : A L<Bio::SeqFeature::SegmentI> object
 Args    : see below
 Status  : Public

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
  my $self = shift;
  my ( $types, $absolute, $baserange, $ranges );
  my %args;
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $types, $absolute, $baserange, $ranges );
      rearrange(
        [ [ qw( TYPE TYPES ) ],
          [ qw( ABSOLUTE ABS ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ]
        ],
        @_
      );
    %args = %_;
  } else {
    ## Types.
    $types = \@_;
    %args = { -types => $types };
  }
  # If -baserange is given but is not a RangeI then the
  # -baserange argument to the superclass should be
  # $self->abs_seq_id().  If -baserange is not given then the argument
  # to the superclass method should be $self unless -absolute is given
  # and true or absolute() is true.
  if( not defined( $baserange ) ) {
    if( $absolute || $self->absolute() ) {
      $baserange = $self->abs_seq_id();
    } else {
      $baserange = $self;
    }
  } elsif( defined( $baserange ) &&
           ( not ref( $baserange ) || not $baserange->isa( 'Bio::RangeI' ) )
         ) {
    $baserange = $self->abs_seq_id();
  }
  $args{ '-baserange' } = $baserange;
  # Oh yeah and make sure it's the only baserange argument..
  $args{ '-base_range' } = $args{ '-baselocation' } =
    $args{ '-base_location' } = $args{ '-base' } =
      $args{ 'baserange' } = $args{ 'base_range' } =
        $args{ 'baselocation' } = $args{ 'base_location' } =
          $args{ 'base' } = undef;

  # Okay so now we just need to make sure that the segment returned is
  # a Segment after all; we do this by overriding the _create_collection
  # method with our own _create_segment method.
  # So we're done!
  return $self->SUPER::get_collection( %args );
} # get_collection(..)

=head2 _create_collection

 Title   : _create_collection
 Usage   : my $segment = $provider->_create_collection(
                              \@args_to_get_collection,
                              @features
                            );
 Function: Factory method for instantiating a collection.
 Args    : a ref to the args used by the get_collection method, and some
           L<Bio::SeqFeatureI> objects to add to the new collection.
 Returns : a new L<Bio::SeqFeature::SegmentI> object
 Status  : Protected

   NOTE THIS CONSTRAINT: Because of our hacky implementation of
get_collection(), we require that the L<Bio::SeqFeature::SegmentI>
that is returned by this method creates segments of its own type
when its get_collection() method is called.  This is the likely
behavior anyway, but it should be noted, just in case.

  This delegates to _create_segment(..)
=cut

# This implementation just delegates to _create_segment.
sub _create_collection {
  shift->_create_segment( @_ );
} # _create_collection

=head2 _create_segment

 Title   : _create_segment
 Usage   : my $segment = $provider->_create_segment(
                              \@args_to_get_collection,
                              @features
                            );
 Function: Factory method for instantiating a segment.
 Args    : a ref to the args used by the get_collection method, and some
           L<Bio::SeqFeatureI> objects to add to the new collection.
 Returns : a new L<Bio::SeqFeature::SegmentI> object
 Status  : Protected

   A SegmentI's get_collection() method must return a SegmentI, and
the returned segment must have its seq_id and start, end, and strand
set appropriately.  The seq_id must be the one used in the
get_collection search; by default it will be this SegmentI's values.
If a range was given then it will be the range of the resulting
SegmentI.  If multiple ranges were given then the union will be used.

   NOTE THIS CONSTRAINT: Because of our hacky implementation of
get_collection(), we require that the L<Bio::SeqFeature::SegmentI>
that is returned by this method creates segments of its own type
when its get_collection() method is called.  This is the likely
behavior anyway, but it should be noted, just in case.

=cut

# This implementation uses the new method of this package, whatever
# this package may be (but it's sure to be some sort of SegmentI..)
sub _create_segment {
  my $self = shift;
  my $args = shift;
  my @features = @_;

  my ( $absolute, $baserange, $ranges );
  if( scalar( @$args ) && $args->[ 0 ] =~ /^-/ ) {
    ( $absolute, $baserange, $ranges );
      rearrange(
        [ [ qw( ABSOLUTE ABS ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ]
        ],
        @_
      );
  }
  # If -baserange is given but is not a RangeI then the
  # -seq_id should be $self->abs_seq_id().  If -baserange is not given
  # then it should be $self unless -absolute is given and true or
  # absolute() is true.
  if( not defined( $baserange ) ) {
    if( $absolute || $self->absolute() ) {
      $baserange = $self->abs_seq_id();
    } else {
      $baserange = $self;
    }
  } elsif( defined( $baserange ) &&
           ( not ref( $baserange ) || not $baserange->isa( 'Bio::RangeI' ) )
         ) {
    $baserange = $self->abs_seq_id();
  }
  my $union_range;
  if( defined( $ranges ) && ref( $ranges ) ) {
    if( ref( $ranges ) eq 'ARRAY' ) {
      $union_range = Bio::RelRange->union( @$ranges );
    } else { # It's gotta be a 'Bio::RangeI'
      $union_range = $ranges;
    }
  }
  my %args_to_new = { '-seq_id' => $baserange };
  if( $baserange ) {
    my $abs_baserange = absSeqId( $baserange );
    if( $abs_baserange && ( $baserange eq $abs_baserange ) ) {
      # $baserange is absolute already..
      $args_to_new{ '-start' } = $union_range->start();
      $args_to_new{ '-end' } = $union_range->end();
      $args_to_new{ '-strand' } = $union_range->strand();
    } else {
      # $baserange is relative, so we need to relativize the union range.
      my $rel_baserange;
      if( $baserange->isa( 'Bio::RelRangeI' ) ) {
        $rel_baserange = $baserange;
      } else {
        $rel_baserange = Bio::RelRange->new( $baserange );
      }
      $args_to_new{ '-start' } =
        $rel_baserange->rel2abs( $union_range->start() );
      $args_to_new{ '-end' } =
        $rel_baserange->rel2abs( $union_range->end() );
      $args_to_new{ '-strand' } =
        $rel_baserange->rel2abs_strand( $union_range->strand() );
    }
  } else {
    # No $baserange.  Oh well.
    $args_to_new{ '-start' } = $union_range->start();
    $args_to_new{ '-end' } = $union_range->end();
    $args_to_new{ '-strand' } = $union_range->strand();
  }
  return $self->new( %args_to_new );
} # _create_segment(..)

# Given a range and a feature, return a feature (maybe the given one)
# that is in coordinates relative to the range.  If the range and the
# feature are on different sequences, an exception will be thrown.
# The range may be a sequence id instead of a Bio::RangeI object.  No
# given feature will be modified at all, so new features will be
# created when a change is necessary.  New features are created using
# the _create_feature_view(..) method, which creates a relativized
# view onto an existing feature (so non-RangeI changes to the returned
# feature will always make a corresponding change to the original one).
sub _relativize_feature {
  my $self = shift;
  my ( $baserange, $feature ) = @_;

  my $abs_baserange = absSeqId( $baserange );
  # Our strategy if the coords are relative to a different range is
  # to get absolute coords and then rerelativize them to the
  # baserange.
  my $abs_feature_range;
  if( $feature->seq_id() ) {
    if( $feature->seq_id() eq $baserange ) {
      return $feature; # It's already ready.
    }
    my $abs_seq_id = $feature->abs_seq_id();
    unless( $abs_seq_id eq $abs_baserange ) { ## Assertion
      $self->throw( "Internal error: a SeqFeatureI object that was to be returned by the features() method of a SimpleSegment is defined over the wrong sequence (it's on '$abs_seq_id' -- we expected it to be on '$abs_baserange').  This is not your fault unless you are the programmer of SimpleSegment or a class on which it is dependent." );
    }
    if( $feature->seq_id() eq $abs_seq_id ) {
      # The feature is already in absolute coords.
      $abs_feature_range = $feature;
    } else {
      $abs_feature_range =
        Bio::RelRange->new(
          '-seq_id' => $feature->seq_id(),
          '-start' => $feature->start(),
          '-end' => $feature->end(),
          '-strand' => $feature->strand()
        );
      # Now absolutify it.
      $abs_feature_range->absolute( 1 );
    }
  } else {
    # If it doesn't have a seq_id, assume it's already absolute.
    if( $abs_baserange eq $baserange ) {
      # If abs_baserange & $baserange are the same, then it's already ready.
      return $feature;
    }
    # The feature is already in absolute coords. 
    $abs_feature_range = $feature;
  } # End if $feature has a $seq_id .. else ..
  # Now it needs to be relativized to $baserange.
  # We can use another range that is relative to $baserange to
  # convert the coords.  Anything will do, so long as it starts at
  # 1 and is on the + strand.
  my $temp_rel_range =
    Bio::RelRange->new(
      '-seq_id' => $baserange,
      '-start' => 1,
      '-end' => 1, # we could use $baserange->length, but it's unnecessary.
      '-strand' => 1
    );
  return $self->_create_feature_view(
    $feature,
    {
      '-seq_id' => $baserange,
      '-start' =>
        $temp_rel_range->abs2rel( $abs_feature_range->start() ),
      '-end' =>
        $temp_rel_range->abs2rel( $abs_feature_range->end() ),
      '-strand' =>
        $temp_rel_range->abs2rel_strand( $abs_feature_range->strand() )
    }
  );
} # _relativize_feature(..)

# Create a feature just like the one given only with a different range
# (the new values are passed in a hash).  If the given feature is
# presently in absolute mode, the new view will be too, but all given
# range values will be interpreted in relative mode.  In the future
# this will be accomplished easily (when the new SeqFeature model is
# ready) by creating a new feature will the same attribute collection
# but a different location.  For now we do a shallow copy of the hash,
# blessed using its own new() method, then change just those parts of
# the range named in the given hash.

sub _create_feature_view {
  my $self = shift;
  my $copy_from = shift;

  my $new_view = $copy_from->new();
  foreach my $key ( keys %$copy_from ) {
    $new_view->{ $key } = $copy_from->{ $key };
  }
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    # The new view should always start out non-absolute while we update stuff
    my $was_absolute =
      ( $new_view->isa( 'Bio::RelRangeI' ) && $new_view->absolute() );
    if( $was_absolute ) {
      $new_view->absolute( 0 );
    }
    my ( $seq_id, $strand, $start, $end, $length ) =
      $self->_rearrange( [ qw( SEQ_ID
                               STRAND
                               START
                               END
                               LENGTH
                             ) ], @_ );
    if( defined( $seq_id ) ) {
      $new_view->seq_id( $seq_id );
    }
    if( defined( $strand ) ) {
      $new_view->strand( $strand );
    }
    if( defined( $start ) ) {
      $new_view->start( $start );
    }
    if( defined( $end ) ) {
      $new_view->end( $end );
    }
    if( defined( $length ) ) {
      $new_view->length( $length );
    }
    if( $was_absolute ) {
      $new_view->absolute( 1 );
    }
  } # End if there's new values to set for the feature view's range.
  return $new_view;
} # _create_feature_view(..)

# Given a range and a callback sub, return a sub that will make sure that all
# features passed to the original sub have been relativized to that range.
sub _create_callback_wrapper {
  my $self = shift;
  my $baserange = shift;
  my $callback = shift;

  return \sub {
    $callback->( $self->_relativize_feature( $baserange, $_[ 0 ] ), $_[ 1 ] );
  };
} # _create_callback_wrapper

## Inner class ##############################################################
#============================================================================
# Bio::SeqFeature::SimpleSegmentIteratorWrapper: A
# Bio::SeqFeature::IteratorI that wraps another to make sure that all
# features that would be returned are relativized to a given range
# beforehand.
#============================================================================
package Bio::SeqFeature::SimpleSegmentIteratorWrapper;
use Bio::Root::Root;
use Bio::SeqFeature::IteratorI;
use vars qw( @ISA );

@ISA = qw( Bio::Root::Root Bio::SeqFeature::IteratorI );

=head2 new

 Title   : new
 Usage   : $iterator = Bio::SeqFeature::SimpleSegmentIteratorWrapper->new(
                         $baserange,
                         $iterator_to_wrap
                       );
 Function: Instantiates a new iterator that relativizes the features
           returned by the given iterator to the given range.
 Returns : a new Bio::SeqFeature::SimpleSegmentIteratorWrapper
 Args    : A L<Bio::RangeI> or sequence id string and an
           L<Bio::SeqFeature::IteratorI> to wrap.
 Status  : Public

=cut

sub new {
  my $class = shift;
  $class = ref( $class ) if ref( $class );
  my ( $baserange, $iterator_to_wrap ) = @_;
  return bless {
		 '_baserange'  => $baserange,
		 '_peer'  => $iterator_to_wrap
	       }, $class;
} # new(..)

=head2 next_feature

 Title   : next_feature
 Usage   : $seq_feature = $iterator->next_feature()
 Function: returns and removes the next feature from the peer iterator,
           relativized.
 Returns : a Bio::SeqFeatureI, or undef if there are no more
 Args    : none
 Status  : Public

=cut

sub next_feature {
  my $self = shift;

  return $self->_relativize_feature(
    $self->{ '_baserange' },
    $self->{ '_peer' }->next_feature()
  );
} # next_feature()

=head2 has_more_features

 Title   : has_more_features
 Usage   : while( $iterator->has_more_features() ) { do something }
 Function: returns true iff there are features in the peer iterator
 Returns : true iff there are more features in the peer iterator
 Args    : none
 Status  : Public

=cut

sub has_more_features {
  my $self = shift;
  return $self->{ '_peer' }->has_more_features();
} # has_more_features()

#============================================================================
## This is the end of SimpleSegmentIteratorWrapper, an inner class of
## SimpleSegment.
#============================================================================
## End Inner class ##########################################################

1;

__END__
