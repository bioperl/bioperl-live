package Bio::SeqFeature::SimpleSegment;

# $Id$
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

use Bio::RelRange qw( &absSeqId );
use Bio::DB::SimpleSegmentProvider;
use Bio::SeqFeature::SimpleCollection;
use Bio::SeqFeature::SegmentI;
@ISA = qw( Bio::RelRange
           Bio::DB::SimpleSegmentProvider
           Bio::SeqFeature::SimpleCollection
           Bio::SeqFeature::SegmentI );

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'
use Bio::SeqFeature::PositionProxy;

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
          : -features (an optional array ref of L<Bio::SeqFeature::SegmentI> objects)
          : -type     (an optional L<Bio::SeqFeature::TypeI> or string)
          : -parent   (an optional L<Bio::DB::SegmentProviderI> parent object)

=cut

sub new {
  my ( $caller, @args ) = @_;

  # Rely on RelRange's new(..) to do most of the work.
  my $self = $caller->SUPER::new( @args );
  $self->_initialize_simple_segment( @args );
  return $self;
} # new(..)

sub _initialize_simple_segment {
  my $self = shift;
  my @args = @_;

  ## TODO: REMOVE
  #warn "_initialize_simple_segment( ".join( ', ', @args )." )";

  return if( $self->{ '_simple_segment_initialized' } );

  my ( $features, $type, $parent, $rest );
  if( scalar( @args ) && ( $args[ 0 ] =~ /^-/ ) ) {
    ( $features, $type, $parent, $rest ) =
      rearrange( [ [ qw( FEATURES SEGMENTS ) ],
                   'TYPE',
                   [ qw( PARENT FACTORY ) ]
                 ], @args );
  }

  ## TODO: REMOVE
  if( $features && ref( $features ) eq 'ARRAY' && @$features ) {
    #warn "_initialize_simple_segment( ".join( ', ', @args )." ): \$features is [ ".join( ', ', @$features )." ]";
  } else {
    #warn "_initialize_simple_segment( ".join( ', ', @args )." ): \$features is $features";
  }

  ## This bizarre -foooo business is to trick the
  ## SimpleSegmentProvider initializer into recognizing that these are
  ## hash arguments (because rearrange(..)'s $rest output has all of
  ## the key names de-dashified, and it looks for a leading dash).
  $self->_initialize_simple_segment_provider( '-foooooo' => 'baaar', %$rest );

  my $added_features = 0;
  if( $features ) {
    # It could be a reference to a list or maybe it is just a single feature..
    if( ref( $features ) eq 'ARRAY' ) {
      foreach my $feature ( @$features ) {
        unless( ref( $feature ) &&
                ( ref( $feature ) ne 'ARRAY' ) &&
                $feature->isa( "Bio::SeqFeature::SegmentI" ) ) {
          $feature = $self->_create_feature( $feature );
        }
        unless( ref( $feature ) && $feature->isa( "Bio::SeqFeature::SegmentI" ) ) {
          $self->throw( "The given feature $feature is not a Bio::SeqFeature::SegmentI.  It is a ".ref( $feature ) );
        }
        unless( $self->_insert_feature( $feature ) ) {
          $self->throw( "duplicate feature: $feature" );
        }
        $added_features++;
      }
    } else {
      unless( ref( $features ) &&
              $features->isa( "Bio::SeqFeature::SegmentI" ) ) {
        $features = $self->_create_feature( $features );
      }
      if( ref( $features ) && $features->isa( 'Bio::SeqFeature::SegmentI' ) ) {
        unless( $self->_insert_feature( $features ) ) {
          $self->throw( "duplicate feature: $features" );
        }
        $added_features++;
      } else {
        $self->throw( "The given feature $features is not a Bio::SeqFeature::SegmentI.  It is a ".ref( $features ) );
      }
    }
  }

  $self->type( $type ) if defined( $type );
  $self->parent_segment_provider( $parent ) if defined( $parent );

  ## TODO: REMOVE.  Testing.
  #unless( $self->start() ) {
  #  $self->start( 1 );
  #  if( $added_features ) {
  #    $self->end( 1 );
  #    $self->adjust_bounds( 1 ); # The argument means yes, shrink if necessary.
  #    ## TODO: REMOVE
  #    if( $self->end() == 1 ) {
  #      #warn "!==! End of $self is 1 after adjust_bounds.  \$features are [ " . join( ', ', @$features ) . " ]";
  #    }
  #    ## TODO: REMOVE
  #    unless( $self->start() ) {
  #      $self->throw( "adjust_bounds caused the start to become 0." );
  #    }
  #  } elsif( ref( $self->seq_id() ) &&
  #           $self->seq_id()->isa( 'Bio::RangeI' ) ) {
  #    ## TODO: REMOVE
  #    #warn "!==! Setting end of $self to " . $self->seq_id()->length() . ".";
  #    # Take on the seq_id's range.
  #    $self->end( $self->seq_id()->length() );
  #  } else {
  #    ## TODO: REMOVE
  #    #warn "!==! Setting end of $self to 1 because its seq_id, " . $self->seq_id() . " is not a RangeI.";
  #    $self->end( 1 );
  #  }
  #}

  $self->{ '_simple_segment_initialized' }++;
  return $self;
} # _initialize_simple_segment(..)

=head2 new_from_segment

 Title   : new_from_segment
 Usage   : my $new_segment =
             Bio::SeqFeature::SimpleSegment->new_from_segment( $copy_from );
 Function: Create a new Bio::SeqFeature::SimpleSegment object by copying
           values from another SimpleSegment object.
 Returns : A new L<Bio::SeqFeature::SimpleSegment> object
 Args    : Another L<Bio::SeqFeature::SimpleSegment> object
 Status  : Protected

  This is a special copy constructor.  It forces the new segment into
  the L<Bio::SeqFeature::SimpleSegment> package, regardless of the
  package that it is called from.  This causes subclass-specfic
  information to be dropped.

  This also does not copy into the new segment the features held in
  the existing segment.  If you would like the new segment to hold the
  same features you must explicitly add them, like so:
    $new_segment->add_features( $copy_from->features() );

  As a special bonus you may also pass an existing hash and it will be the
  blessed an anointed object that is returned, like so:
    $new_segment =
      Bio::SeqFeature::SimpleSegment->new_from_segment(
        $copy_from,
        $new_segment
      );

=cut

sub new_from_segment {
  my $pack = shift; # ignored
  my $segment = shift || $pack;
  my $new_segment = shift;
  $new_segment = Bio::RelRange->new_from_relrange( $segment, $new_segment );
  $new_segment =
    Bio::SeqFeature::SimpleCollection->new_from_collection(
      $segment,
      $new_segment
    );
  return bless $new_segment, __PACKAGE__;
} # new_from_segment(..)

=head2 _create_feature

 Title   : create_feature
 Usage   : my $new_seq_feature = $segment->_create_feature( $feature_data );
 Function: Factory method for instantiating a SeqFeatureI object.
 Returns : A new L<Bio::SeqFeature::SegmentI> object, or $feature_data.
 Args    : A single argument of any data type (see below).
 Status  : Protected

 The single argument may be of any type.  _create_feature(..) will be
 called by _insert_feature whenever the feature argument to that
 method is not a SeqFeatureI object.  If _create_feature(..) is unable
 to make a feature from the given data it must return the argument it
 was given.

=cut

## This one makes features when given [ $start, $end ] pairs.  If this
## SimpleSegment implements the SeqFeatureI interface then
## $self->new(..) will be used instead of
## Bio::SeqFeature::Generic->new(..).  New features will have the same
## type as $self does and will have $self as their $seq_id.  The given
## start and end positions are interpreted relative to $self.
sub _create_feature {
  my $self = shift;
  my $feature_data = shift;
  if( ref( $feature_data ) eq 'ARRAY' ) {
    # This is for feature data like [ [ start0, end0 ], [ start1, end1 ] ].
    my ( $start, $end ) = @{ $feature_data };
    
    # The following line should be unnecessary.
    #next unless defined $start && defined $end;
    
    # Strandedness defaults to our own, but start > end forces -1.
    my $strand = $self->strand();
    if( $start > $end ) {
      ( $start, $end ) = ( $end, $start );
      $strand = -1;
    }
    if( $self->isa( 'Bio::SeqFeature::SegmentI' ) ) {
      return $self->new( -start  => $start,
                         -end    => $end,
                         -strand => $strand,
                         -type   => $self->type(),
                         -seq_id => $self,
                         -parent => $self ) || $feature_data;
    } else {
      return Bio::SeqFeature::Generic->new(
        -start  => $start,
        -end    => $end,
        -strand => $strand,
        -type   => $self->type(),
        -seq_id => $self,
        -parent => $self
      ) || $feature_data;
    }
  }
  # If we're at this point then we've been unable to make a new feature.
  return $feature_data;
} # _create_feature(..)

=head2 unique_id

 Title   : unique_id
 Usage   : my $unique_id = $segment->unique_id( [$new_unique_id] )
 Function: This is a unique identifier that identifies this segment object.
           If not set, will return undef per L<Bio::LocallyIdentifiableI>
           If a value is given, the unique_id will be set to it, unless that
           value is the string 'undef', in which case the unique_id will
           become undefined.
 Returns : The current (or former, if used as a set method) value of unique_id
 Args    : [optional] a new string unique_id or 'undef'

=cut

sub unique_id {
  my ( $self, $value ) = @_;
  my $current_value = $self->{ '_unique_id' };
  if ( defined $value ) {
    if( !$value || ( $value eq 'undef' ) ) {
      $self->{ '_unique_id' } = undef;
    } else {
      $self->{ '_unique_id' } = $value;
    }
  }
  return $current_value;
} # unique_id()

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

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
    %args = @_;
  } else {
    ## Types.
    $types = \@_;
    %args = ( '-types' => $types );
  }
  # If -baserange is given but is not a RangeI then the
  # -baserange argument to the superclass should be
  # $self->abs_seq_id().  If -baserange is not given then the argument
  # to the superclass method should be $self unless -absolute is given
  # and true or absolute() is true.
  if( not defined( $baserange ) ) {
    if( $absolute || $self->absolute() ) {
      $baserange = $self->abs_range();
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
  ## TODO: What if they're uppercase?
  delete @args{ qw( -base_range -baselocation -base_location -base
                    baserange base_range baselocation base_location base ) };

  # We do pass -absolute up to SimpleCollection.  Although it does not
  # document it, passing -absolute up to SimpleCollection affects the
  # sort order of the returned features. (sorts by abs positions)

  if( $absolute ||
      ( defined( $baserange ) && ref( $baserange ) &&
        $baserange->isa( 'Bio::RelRangeI' ) && $baserange->absolute() ) ) {
    $args{ '-absolute' } = 1;
  }
  # Oh yeah and make sure it's the only absolute argument..
  ## TODO: What if they're uppercase?
  delete @args{ qw( -abs absolute abs ) };

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

  # We also have promised that the returned features will be relative
  # to whatever baserange ends up being.
  if( $args{ '-absolute' } ) {
    # In this case it doesn't matter.
    return $self->SUPER::features( %args );
  }
  # There's 3 options for this: by default features() returns a
  # list. if -iterator is given then it returns an iterator.  if -callback
  # is given then some method is called.  We need to deal with each option
  # separately.
  if( $iterator ) {
    # Wrap the returned iterator in our own..
    my $super_iterator = $self->SUPER::features( %args );
    ## SimpleSegment_IteratorWrapper is defined in this file, below.
    return Bio::SeqFeature::SimpleSegment_IteratorWrapper->new(
      $baserange,
      $super_iterator
    );
  } elsif( $callback ) {
    # Wrap the callback in our own..
    my $new_callback =
      $self->_create_callback_wrapper( $baserange, $callback );
    $args{ '-callback' } = $new_callback;
    # Oh and make sure it's the only callback argument..
    ## TODO: What if they're uppercase?
    delete @args{ qw( -call_back callback call_back ) };
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
  my ( $types, $unique_ids, $names, $attributes, $baserange, $ranges,
       $absolute );
  my %args;
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $types, $unique_ids, $names, $attributes, $baserange, $ranges,
      $absolute ) =
      rearrange(
        [ [ qw( TYPE TYPES ) ],
          [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
          [ qw( NAME NAMES DISPLAY_NAME DISPLAY_NAMES DISPLAYNAME DISPLAYNAMES ) ],
          [ qw( ATTRIBUTE ATTRIBUTES ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ],
          [ qw( ABSOLUTE ABS ) ]
        ],
        @_
      );
    %args = @_;
  } elsif( scalar( @_ ) ) {
    ## Types.
    $types = \@_;
    %args = ( '-types' => $types );
  }

  ## Fix up types.
  if( $types ) {
    unless( ref $types eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $types = [ $types ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$types ) ||
       ( ( scalar( @$types ) == 1 ) && !( $types->[ 0 ] ) )
      ) {
      undef $types;
    }
  }

  ## Fix up unique_ids.
  if( $unique_ids ) {
    unless( ref $unique_ids eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $unique_ids = [ $unique_ids ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$unique_ids ) ||
       ( ( scalar( @$unique_ids ) == 1 ) && !( $unique_ids->[ 0 ] ) )
      ) {
      undef $unique_ids;
    }
  }

  ## Fix up names.
  if( $names ) {
    unless( ref $names eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $names = [ $names ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$names ) ||
       ( ( scalar( @$names ) == 1 ) && !( $names->[ 0 ] ) )
      ) {
      undef $names;
    }
  }

  ## Attributes better be a hash ref if it is anything.
  if( $attributes ) {
    unless( ref( $attributes ) eq 'HASH' ) {
      $self->throw( "The -attributes argument must be a HASH REF." );
    }
  }

  ## Fix up ranges.
  if( $ranges ) {
    unless( ref $ranges eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $ranges = [ $ranges ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$ranges ) ||
       ( ( scalar( @$ranges ) == 1 ) && !( $ranges->[ 0 ] ) )
      ) {
      undef $ranges;
    }
  }

  ## TODO: REMOVE
  #warn "SimpleSegment->get_collection( ".join( ', ', ( my @foo = %args ) ). " )" if Bio::Graphics::Browser::DEBUG;

  # If -baserange is given but is not a RangeI then the
  # -baserange argument to the superclass should be
  # $self->abs_seq_id().  If -baserange is not given then the argument
  # to the superclass method should be $self unless -absolute is given
  # and true or absolute() is true.
  if( not defined( $baserange ) ) {
    if( $absolute || $self->absolute() ) {
      $baserange = $self->abs_range();
      $absolute = 1;
    } else {
      $baserange = $self;
    }
  } elsif( defined( $baserange ) &&
           ( not ref( $baserange ) || not $baserange->isa( 'Bio::RangeI' ) )
         ) {
    $baserange = $self->abs_range();
  }
  $args{ '-baserange' } = $baserange;
  # Oh yeah and make sure it's the only baserange argument..
  ## TODO: What if they're uppercase?
  delete @args{ qw( -base_range -baselocation -base_location -base
                    baserange base_range baselocation base_location base ) };

  # Okay so now we just need to make sure that the segment returned is
  # a Segment after all; we do this by overriding the _create_collection
  # method with our own _create_segment method.
  # So we're done!
  ## TODO: REMOVE
  #warn "SimpleSegment->Bio::SeqFeature::SimpleCollectionProvider::get_collection( ".join( ', ', ( my @foo = %args ) ). " )" if Bio::Graphics::Browser::DEBUG;

  my $segment =
    $self->Bio::SeqFeature::SimpleCollectionProvider::get_collection( %args );

  ## TODO: REMOVE. Testing.
  if( defined( $segment ) && !ref( $segment ) ) {
    $self->throw( "Geez.  The result of SimpleSegment->Bio::SeqFeature::SimpleCollectionProvider::get_collection( ".join( ', ', ( my @foo = %args ) ). " ) is $segment, a ".ref( $segment )."." );
  }

  # If there's no features, AND the request named something in
  # particular, then we have to return undef to indicate that the
  # request failed.
  if(
     ( !defined( $segment ) || ( $segment->feature_count() == 0 ) ) &&
     ( $types || $unique_ids || $names || $attributes || $ranges )
    ) {
    return ( wantarray ? () : undef );
  } else {
    return $segment;
  }
} # get_collection(..)

=head2 _create_collection

 Title   : _create_collection
 Usage   : my $segment = $provider->_create_collection(
                           \@args_to_get_collection,
                           @features
                         );
           or
           my $segment = $provider->_create_collection(
                           \@args_to_get_collection,
                           'lookup'
                         );
 Function: Factory method for instantiating a segment.
 Args    : a ref to the args used by the get_collection method, and some
           L<Bio::SeqFeature::SegmentI> objects to add to the new
           segment, or 'lookup', meaning that the args should be used
           to add the features to the new collection.
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
  return shift->_create_segment( @_ );
} # _create_collection

=head2 _create_segment

 Title   : _create_segment
 Usage   : my $segment = $provider->_create_segment(
                           \@args_to_get_collection,
                           @features
                         );
           or
           my $segment = $provider->_create_segment(
                           \@args_to_get_collection,
                           'lookup'
                         );
 Function: Factory method for instantiating a segment.
 Args    : a ref to the args used by the get_collection method, and some
           L<Bio::SeqFeature::SegmentI> objects to add to the new
           segment, or 'lookup', meaning that the args should be used
           to add the features to the new collection.
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
  ## HACK because sometimes the args are passed in as a hash and
  ## sometimes as a list.
  if( $args && ( ref( $args ) eq 'HASH' ) ) {
    my @args_list = %$args;
    $args = \@args_list;
  }
  my @features = @_;
  if( @features && ( $features[ 0 ] eq 'lookup' ) ) {
    @features = $self->features( @$args );
  }

  ## TODO: REMOVE
  #warn "SimpleSegment::_create_segment( { ".join( ', ', @$args )." }, ( ".join( ', ', @features )." )" if Bio::Graphics::Browser::DEBUG;

  ## TODO: REMOVE
  if( @features ) {
    #warn "creating segment with these features: ( ", join( ', ', @features ), " )" if Bio::Graphics::Browser::DEBUG;
    #warn "The first one of those has " . $features[ 0 ]->feature_count() . " features." if Bio::Graphics::Browser::DEBUG;
  } else {
    #warn "creating segment with no features.\n";
  }

  my ( $absolute, $baserange, $ranges, $carryover_args );
  if( scalar( @$args ) && $args->[ 0 ] =~ /^-/ ) {
    ( $absolute, $baserange, $ranges, $carryover_args ) =
      rearrange(
        [ [ qw( ABSOLUTE ABS ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ]
        ],
        @$args
      );
    if( $carryover_args ) {
      # Unfortunately it changes the rest of the args, so let's make
      # them nice again.
      my @ca_list = %$carryover_args;
      for( my $i = 0; $i < $#ca_list; $i += 2 ) {
        $ca_list[ $i ] = '-'.lc( $ca_list[ $i ] );
      }
      %$carryover_args = @ca_list;
    }
  }
  # If -baserange is given but is not a RangeI then the
  # -seq_id should be $self->abs_seq_id().  If -baserange is not given
  # then it should be $self unless -absolute is given and true or
  # absolute() is true.
  if( !defined( $baserange ) ) {
    if( $absolute || $self->absolute() ) {
      $baserange = $self->abs_range();
    } else {
      $baserange = $self;
    }
  } elsif( defined( $baserange ) &&
           ( !ref( $baserange ) || !$baserange->isa( 'Bio::RangeI' ) )
         ) {
    ## TODO: REMOVE?
    warn "The baserange that was given, $baserange, isa ".ref( $baserange ).", but we need a RangeI.  Using \$self->abs_range().";
    $baserange = $self->abs_range();
  }
  my $union_range;
  if( defined( $ranges ) && ref( $ranges ) ) {
    if( ref( $ranges ) eq 'ARRAY' ) {
      $union_range = Bio::RelRange->union( @$ranges );
    } else { # It's gotta be a 'Bio::RangeI'
      $union_range = $ranges;
    }
  } elsif( @features ) {
    $union_range = Bio::RelRange->union( @features );
  } else {
    $union_range = $self;
  }
  ## TODO: REMOVE
  #warn "SimpleSegment::_create_segment(..): union_range is $union_range, baserange is " . $baserange->Bio::RelRangeI::toString() . '.' if Bio::Graphics::Browser::DEBUG;
  my %args_to_new =
    (
     '-seq_id'   => $baserange,
     '-parent'   => $self,
     '-features' => \@features,
     ( ( defined $carryover_args ) ? ( my @foo = %$carryover_args ) : () )
    );
  if( $baserange ) {
    ## TODO: Put back
    #my $abs_baserange = $baserange->abs_seq_id();
    ## TODO: REMOVE.  Testing.
    my $abs_baserange = eval{ $baserange->abs_range() } || $baserange;
    if( $abs_baserange && ( $baserange eq $abs_baserange ) ) {
      ## TODO: REMOVE
      #warn "baserange $baserange is absolute already.";
      # $baserange is absolute already..
      $args_to_new{ '-start' }  = $union_range->start();
      $args_to_new{ '-end' }    = $union_range->end();
      $args_to_new{ '-strand' } = $union_range->strand();
    } else {
      ## TODO: REMOVE
      #warn "baserange $baserange is not absolute.  \$abs_baserange is $abs_baserange.";
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
    $args_to_new{ '-start' }  = $union_range->start();
    $args_to_new{ '-end' }    = $union_range->end();
    $args_to_new{ '-strand' } = $union_range->strand();
  }
  ## TODO: REMOVE
  #warn "SimpleSegment::_create_segment(..): \%args_to_new are { ".join( ', ', ( my @foo = %args_to_new ) )." }" if Bio::Graphics::Browser::DEBUG;
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

  ## TODO: REMOVE
  #warn "SimpleSegment::_relativize_feature( $baserange, $feature )" if Bio::Graphics::Browser::DEBUG;

  my $abs_baserange = absSeqId( $baserange );
  # Our strategy if the coords are relative to a different range is
  # to get absolute coords and then rerelativize them to the
  # baserange.
  my $abs_feature_range;
  if( $feature->seq_id() ) {
    if( $feature->seq_id() eq $baserange ) {
      ## TODO: REMOVE
      #warn "                                                $feature" if Bio::Graphics::Browser::DEBUG;
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
      $abs_feature_range = Bio::RelRange->new( $feature );
      # Now absolutify it.
      $abs_feature_range->absolute( 1 );
    }
  } else {
    # If it doesn't have a seq_id, assume it's already absolute.
    if( $abs_baserange eq $baserange ) {
      # If abs_baserange & $baserange are the same, then it's already ready.
      ## TODO: REMOVE
      #warn "                                                $feature" if Bio::Graphics::Browser::DEBUG;
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

  my $new_start = $temp_rel_range->abs2rel( $abs_feature_range->start() );
  my $new_end = $temp_rel_range->abs2rel( $abs_feature_range->end() );
  my $new_strand =
    $temp_rel_range->abs2rel_strand( $abs_feature_range->strand() );
  # Don't let the strand get forced to negative:
  if( ( $new_strand >= 0 ) && ( $new_end < $new_start ) ) {
    ( $new_start, $new_end ) = ( $new_end, $new_start );
  }
  #return $self->_create_feature_view(
  #  $feature,
  #  (
  #    '-seq_id' => $baserange,
  #    '-start' => $new_start,
  #    '-end' => $new_end,
  #    '-strand' => $new_strand
  #  )
  #);
  my $new_feature = $self->_create_feature_view(
    $feature,
    (
      '-seq_id' => $baserange,
      '-start' => $new_start,
      '-end' => $new_end,
      '-strand' => $new_strand
    )
  );
  ## TODO: REMOVE
  #warn "                                                $new_feature" if Bio::Graphics::Browser::DEBUG;
  return $new_feature;
} # _relativize_feature(..)

# Create a feature just like the one given only with a different range
# (the new values are passed in a hash).  If the given feature is
# presently in absolute mode, the new view will be too, but all given
# range values will be interpreted in relative mode.  In the future
# this will be accomplished easily (when the new SeqFeature model is
# ready) by creating a new feature will the same attribute collection
# but a different location.  For now we use a PositionProxy.
sub _create_feature_view {
  my $self = shift;
  my $copy_from = shift;

  return
    Bio::SeqFeature::PositionProxy->new(
      '-peer' => $copy_from,
      '-absolute' => $copy_from->absolute(),
      @_
    );
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
# Bio::SeqFeature::SimpleSegment_IteratorWrapper: A
# Bio::SeqFeature::IteratorI that wraps another to make sure that all
# features that would be returned are relativized to a given range
# beforehand.
#============================================================================
package Bio::SeqFeature::SimpleSegment_IteratorWrapper;
use Bio::Root::Root;
use Bio::SeqFeature::IteratorI;
use vars qw( @ISA );

@ISA = qw( Bio::Root::Root Bio::SeqFeature::IteratorI );

=head2 new

 Title   : new
 Usage   : $iterator = Bio::SeqFeature::SimpleSegment_IteratorWrapper->new(
                         $baserange,
                         $iterator_to_wrap
                       );
 Function: Instantiates a new iterator that relativizes the features
           returned by the given iterator to the given range.
 Returns : a new Bio::SeqFeature::SimpleSegment_IteratorWrapper
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
 Returns : a Bio::SeqFeature::SegmentI, or undef if there are no more
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
## This is the end of SimpleSegment_IteratorWrapper, an inner class of
## SimpleSegment.
#============================================================================
## End Inner class ##########################################################

1;

__END__
