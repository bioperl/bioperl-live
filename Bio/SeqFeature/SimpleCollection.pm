package Bio::SeqFeature::SimpleCollection;

# $Id$
## A simple implementation of the Bio::SeqFeature::CollectionI interface.

=head1 NAME

Bio::SeqFeature::SimpleCollection - An in-memory implementation of the
CollectionI interface.

=head1 SYNOPSIS

 use Bio::SeqFeature::SimpleCollection;
 my $collection = new Bio::SeqFeature::SimpleCollection;
 $collection->add_features( @featurelist );
 
 $collection->features(-location => new Bio::Location::Simple
                                        (-start=> 1, -end => 300),
                       -rangetype=>'overlaps');

=head1 DESCRIPTION

A SimpleCollection is an in-memory store of SeqFeatureIs that
implements the CollectionI interface.  It supports updating, adding,
and removing SeqFeatures in its store, although explicit updating is
unnecessary as the SeqFeatureIs returned by the get_collection and
features methods are the same as those stored herein (ie. there is no
external backing store).

Features can be filtered by the following attributes:

  1) their location, perhaps relative to a range (with a choice
     between overlapping, contained within, or completely containing a
     range).

  2) their type

  3) other attributes using tag/value semantics

Access to the contained features can be achieved using the
CollectionProviderI interface to produce a new semantic view on this
collection, or via accessing the feature list directly.  Access to the
feature list can be achieved using multiple techniques:

  1) as another Collection

  2) fetching entire list of features at a time

  3) fetching an iterator across features

  4) a callback

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
use vars qw( $VERSION @ISA );
use overload 
  '""' => 'toString',
  cmp   => '_cmp';

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'

require Bio::SeqFeature::Generic;
require Bio::SeqFeature::SimpleIterator;
use Bio::SeqFeature::SimpleCollectionProvider;
use Bio::SeqFeature::CollectionI;
use Bio::Location::Simple;

$VERSION = '0.01';
@ISA = qw( Bio::SeqFeature::SimpleCollectionProvider
           Bio::SeqFeature::CollectionI );

=head2 new

 Title   : new
 Usage   : my $obj =
             new Bio::SeqFeature::SimpleCollection( @features );
 Function: Builds a new Bio::SeqFeature::SimpleCollection object 
 Returns : Bio::SeqFeature::SimpleCollection
 Args    : SeqFeatureI objects to store herein

=cut

sub new {
  my( $class, @args ) = @_;

  my $self = $class->SUPER::new( @args );

  return $self;
} # new(..)

=head2 new_from_collection

 Title   : new_from_collection
 Usage   : my $new_collection =
             Bio::SeqFeature::SimpleCollection->new_from_collection(
               $copy_from
             );
 Function: Create a new Bio::SeqFeature::SimpleCollection object by copying
           values from another SimpleCollection object.
 Returns : A new L<Bio::SeqFeature::SimpleCollection> object
 Args    : Another L<Bio::SeqFeature::SimpleCollection> object
 Status  : Protected

  This is a special copy constructor.  It forces the new collection into
  the L<Bio::SeqFeature::SimpleCollection> package, regardless of the
  package that it is called from.  This causes subclass-specfic
  information to be dropped.

  This also does not copy into the new collection the features held in
  the existing collection.  If you would like the new collection to hold the
  same features you must explicitly add them, like so:
    $new_collection->add_features( $copy_from->features() );

  As a special bonus you may also pass an existing hash and it will be the
  blessed an anointed object that is returned, like so:
    $new_collection =
      Bio::SeqFeature::SimpleCollection->new_from_collection(
        $copy_from,
        $new_collection
      );

=cut

sub new_from_collection {
  my $pack = shift; # ignored
  my $collection = shift || $pack;
  my $new_collection = shift;
  $new_collection =
    Bio::SeqFeature::SimpleCollectionProvider->new_from_collectionprovider(
      $collection,
      $new_collection
    );
  @{ $new_collection }{ qw( _sorted ) } =
    @{ $collection }{ qw( _sorted ) };
  return bless $new_collection, __PACKAGE__;
} # new_from_collection(..)

=head2 sorted

 Title   : sorted
 Usage   : my $is_sorted = $collection->sorted( [$new_sorted] );
 Function: Getter/setter for the sorted flag.  When this flag is true the
           features() method always returns features in sorted order.
 Returns : The current (or former, if used as a setter) value of the sorted
           flag.
 Args    : (optional) a new value for the sorted flag
 Status  : Public

  See the features() method for more information on sorting.

  Please note that setting this flag may cause the internal
representation of the features to change, so if you eval{ $c->sorted(
1 ); $c->sorted( 0 ); }, the order afterwards is not guaranteed to be
what it was beforehand.  The SimpleCollection implementation does not
cause an internal representation change, but subclasses might.

=cut

sub sorted {
  my ( $self, $new_sorted ) = @_;
  my $current_val = $self->{ '_sorted' };
  if( defined $new_sorted ) {
    $self->{ '_sorted' } = $new_sorted;
  }
  return $current_val;
} # sorted()

=head2 add_features

 Title   : add_features
 Usage   : my @added = $collection->add_features( @feature_list );
 Function: Adds the given features to this Collection.
 Returns : The features added (or their count, in scalar context).
 Args    : An array of L<Bio::SeqFeature::SegmentI>s
 Status  : Public

=cut

sub add_features {
  my $self = shift;
  my @feature_list = @_;

  my $count = 0;
  my @added;
  foreach my $feature ( @feature_list ) {
    if( $feature = $self->_insert_feature( $feature ) ) {
      if( wantarray ) {
        push( @added, $feature );
      } else {
        $count++;
      }
    }
  }

  return ( wantarray ? @added : $count );
} # add_features(..)

=head2 remove_features

 Title   : remove_features
 Usage   : my @removed = $collection->remove_features( @feature_list )
 Function: Removes the requested sequence features
 Returns : The removed features (or their count, in scalar context)
 Args    : An array of L<Bio::SeqFeature::SegmentI>s or their unique_ids
 Status  : Public

=cut

sub remove_features {
  my $self = shift;
  my @feature_list = @_;

  my $count = 0;
  my @removed;
  foreach my $feature ( @feature_list ) {
    ## Special case, to allow the list to include unique_ids.
    if( ref \$feature eq 'STRING' ) {
      $feature = Bio::SeqFeature::Generic->new( '-unique_id' => $feature );
    }
    if( $feature = $self->_remove_feature( $feature ) ) {
      if( wantarray ) {
        push( @removed, $feature );
      } else {
        $count++;
      }
    }
  }
  return ( wantarray ? @removed : $count );
} # remove_features()

=head2 features

 Title   : features
 Usage   : @features = $collection->features( %args );
           OR
           @features = $collection->features( @types );
 Returns : a list of L<Bio::SeqFeature::SegmentI> objects,
           OR
           (when the -iterator option is true) an L<Bio::SeqFeature::IteratorI>
           OR
           (when the -callback argument is given) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public

This routine will retrieve features associated with this collection
object.  It can be used to return all features, or a subset based on
their type, location, or attributes.

If ranges are specified using the -ranges argument, then these ranges
will be used to narrow the results, according to the specified
-rangetype and -strandtype arguments.

If a -baserange is specified or is provided by default* then
unqualified ranges given with the -ranges argument will be interpreted
as relative to that baserange.  Note that this only applies to
unqualified ranges, ie. ranges that have no defined seq_id.  You may
force absolute range interpretations by giving a -baserange that is
not a L<Bio::RangeI> (such as the string 'absolute') or by qualifying all
given ranges.

Footnote (*): All implementing classes that also implement L<Bio::RangeI>
              B<must> provide $self as the default baserange!

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
                 the -range argument is relative.  There may be a
                 default -baserange.  If this CollectionI is also a
                 L<Bio::RangeI>, then the default -baserange should be
                 itself.  Note that the baserange affects the sort order.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch, -rangetype, and -baserange.
  -ranges        An array reference to multiple ranges.

  -rangetype     One of "overlaps", "contains", or "contained_in".

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
arguments consisting of the current feature and this CollectionI
object, and must return a true value. If the code returns a false
value, feature retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

-callback and -sort are mutually exclusive options.  If -sort is
defined, then -callback is ignored.  If you want to do a sorted
callback, set the sorted() flag of this collection to true.

If -sort or sorted() is true then the features will be returned in
order of the features' start positions.  This order will be reversed
if the baserange has a negative strand (remember that a CollectionI
implementation that is also a L<Bio::RangeI> must provide itself as
the default baserange, but this may be overridden by the -baserange
argument).

Note that no guarantees are made by the CollectionI interface about
the order of the features, except when the sorted() flag is true or
when the -sort option is given to the features method.  Therefore
the implementation may choose to reorder the underlying data structure
to better accomodate -sorted feature requests as a result of a
features() call.  When this happens the CollectionI's sorted() flag
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

  my ( $types, $unique_ids, $namespace, $names, $attributes, $baserange, $ranges, $strandmatch, $rangetype, $iterator, $callback, $sort, $absolute );
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $types, $unique_ids, $namespace, $names, $attributes, $baserange, $ranges, $strandmatch, $rangetype, $iterator, $callback, $sort, $absolute ) =
      rearrange(
        [ [ qw( TYPE TYPES ) ],
          [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
          [ qw( NAMESPACE NAME_SPACE CLASS ) ],
          [ qw( NAME NAMES DISPLAY_NAME DISPLAY_NAMES DISPLAYNAME DISPLAYNAMES ) ],
          [ qw( ATTRIBUTE ATTRIBUTES ) ],
          [ qw( BASERANGE BASE_RANGE BASELOCATION BASE_LOCATION BASE ) ],
          [ qw( RANGE RANGES LOCATION LOCATIONS LOC ) ],
          [ qw( STRANDMATCH STRAND_MATCH ) ],
          [ qw( RANGETYPE RANGE_TYPE ) ],
          [ qw( ITERATOR STREAM ) ],
          [ qw( CALLBACK CALL_BACK ) ],
          [ qw( SORT SORTED ) ],
          [ qw( ABSOLUTE ABS ) ]
        ],
        @_
      );
  } else {
    ## Types.
    $types = \@_;
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

  ## We can use ourselves as the baserange, if we are a Bio::RangeI.
  if( !defined( $baserange ) && $self->isa( 'Bio::RangeI' ) ) {
    $baserange = $self;
  }

  ## Derelativize, man.
  if( $ranges && @$ranges && $baserange ) {
    my @new_ranges;
    foreach my $range ( @$ranges ) {
      unless( ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
        $self->throw( "Expected the -ranges argument to be a reference to a list of Bio::RangeI objects, but it contains something incompatible: " . overload::StrVal( $range ) );
      }
      unless( defined $range->seq_id() ) {
        $range = Bio::RelRange->new(
                   -seq_id => $baserange,
                   -start =>  $range->start(),
                   -end =>    $range->end(),
                   -strand => $range->strand()
                 );
      } elsif( !$range->isa( 'Bio::RelRangeI' ) ) {
        $range = Bio::RelRange->new(
                   -seq_id => $range->seq_id(),
                   -start =>  $range->start(),
                   -end =>    $range->end(),
                   -strand => $range->strand()
                 );
      }
      $range->absolute( 1 );
      push( @new_ranges, $range );
    }
    $ranges = \@new_ranges;
  } # End derelativizing ranges

  ## -strandmatch defaults to 'weak'
  unless( defined $strandmatch ) {
    $strandmatch = 'weak';
  }

  ## -rangetype defaults to 'overlaps'
  unless( defined $rangetype ) {
    $rangetype = 'overlaps';
  }

  ## -iterator and -callback are mutually exclusive.
  if( $iterator && $callback ) {
    $self->throw( "The -iterator and -callback options are mutually exclusive, and yet both have been specified.  Perhaps you could apply your callback method to each element returned by the iterator?" );
  }

  ## -sort is redundant if sorted() is true.
  if( $sort && $self->sorted() ) {
    undef $sort;
  }

  ## -sort and -callback are mutually exclusive.
  if( $sort && $callback ) {
    $self->throw( "The -sort and -callback options are mutually exclusive, and yet both have been specified.  Perhaps you could first set the sorted() property of this collection to true, and then try again." );
  }

  ## We actually implement the sorting on-the-fly anyway, so the above
  ## redundancy check is silly (but we keep it for the mutex check)
  $sort ||= $self->sorted();
  ## Do a reverse sort iff the baserange has a negative strand.
  my $reverse_sort = 0;
  if( $baserange && ref( $baserange ) && $baserange->isa( 'Bio::RangeI' ) ) {
    if( $absolute && $baserange->isa( 'Bio::RelRangeI' ) ) {
      $reverse_sort = ( $baserange->abs_strand() < 0 );
    } elsif( !$absolute ) {
      $reverse_sort = ( $baserange->strand() < 0 );
    }
  }
  ## This may be horribly inefficient, but what can ya do?  That's why
  ## this is called SimpleCollection.  Make your own if you want
  ## something better.
  my @features_to_return;

  ## This next line is fancy speak for "for each feature we've got,
  ## perhaps sorted, perhaps reverse-sorted.."
  foreach my $feature 
    (
     ( $sort ?
       sub{
         sort { ( $reverse_sort ?
                  ( $absolute ?
                    $b->abs_high() <=> $a->abs_high() :
                    $b->high() <=> $a->high() ) :
                  ( $absolute ?
                    $a->abs_low() <=> $b->abs_low() :
                    $a->low() <=> $b->low() ) ) } @_;
       } :
       sub { return @_; }
     )->( map { ( values %{ $self->{ '_seq_id_to_feature_table' }->
                                   { $_ }->
                                   { '_identifiable_features' } } ),
                  ( map @{ $_ },
                    ( values %{ $self->{ '_seq_id_to_feature_table' }->
                                       { $_ }->
                                       { '_anonymous_features' } } ) ) }
               keys %{ $self->{ '_seq_id_to_feature_table' } } )
    ) {

    ## TODO: REMOVE
    unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeature::SegmentI' ) ) {
      print STDERR "--------------Hey, man.  $feature ain't a feature!=======\n";
    }

    ## Filter them:
    my $passes_filter = 1;

    ## Filter on types
    if( $passes_filter && $types ) {
      # Failure until proven success:
      $passes_filter = 0;
      my $feature_type = $feature->type();
      foreach my $type ( @$types ) {
        if( ref( $feature_type ) &&
            $feature_type->isa( 'Bio::DB::GFF::Typename' ) &&
            $feature_type->match( $type ) ) {
          # success
          $passes_filter = 1;
          last;
        } elsif( ref( $type ) && $type->isa( 'Bio::SeqFeature::TypeI' ) ) {
          if( ( $type eq $feature_type ) ||
              $type->is_descendent( $feature_type ) ) {
            # success
            $passes_filter = 1;
            last;
          }
        } elsif( ref( $type ) && $type->isa( 'Bio::DB::GFF::Typename' ) ) {
          ## I've added this to support the Typename stuff used by GFF.
          if( $type->match( $feature_type ) ) {
            # success
            $passes_filter = 1;
            last;
          }
        } else {
          ## $type is a string.
          if( ref( $feature_type ) &&
              $feature_type->isa( 'Bio::SeqFeature::TypeI' ) ) {
            if( ( $feature_type eq $type ) ||
                $feature_type->is_ancestor( $type ) ) {
              # success
              $passes_filter = 1;
              last;
            }
          } else {
            ## $feature_type is also a string.
            if( $type =~ /$feature_type/i ) { # ignore case (should we? ok)
              # success
              $passes_filter = 1;
              last;
            }
          }
        } # End if $type is a TypeI .. else ..
      } # End for each $type
    } # End if we're filtering on types

    ## Filter on unique_ids
    ## If they've provided a $namespace, also try "$namespace:$feature_id"
    ## and "$namespace:$id"
    if( $passes_filter && $unique_ids ) {
      # Failure until proven success:
      $passes_filter = 0;
      my $feature_id = $feature->unique_id() || $feature->display_name();
      foreach my $unique_id ( @$unique_ids ) {
        if( ( $feature_id eq $unique_id ) ||
            ( $namespace && ( ( $namespace.':'.$feature_id ) eq $unique_id ) ) ||
            ( $namespace && ( $feature_id eq ( $namespace.':'.$unique_id ) ) )
          ) {
          # success
          $passes_filter = 1;
          last;
        }
      } # End for each $unique_id
    } # End if we're filtering on unique_ids

    ## Filter on names
    ## If they've provided a $namespace, also try "$namespace:$feature_name"
    ## and "$namespace:$name" and "$namespace:$feature_id"
    if( $passes_filter && $names ) {
      # Failure until proven success:
      $passes_filter = 0;
      my $feature_name = $feature->display_name();
      my $feature_id = $feature->unique_id();
      foreach my $name ( @$names ) {
        if( ( $feature_name eq $name ) ||
            ( $namespace && ( ( $namespace.':'.$feature_name ) eq $name ) ) ||
            ( $namespace && ( $feature_name eq ( $namespace.':'.$name ) ) ) ||
            ( $feature_id eq $name ) ||
            ( $namespace && ( ( $namespace.':'.$feature_id ) eq $name ) ) ||
            ( $namespace && ( $feature_id eq ( $namespace.':'.$name ) ) )
          ) {
          # success
          $passes_filter = 1;
          last;
        }
      } # End for each $name
    } # End if we're filtering on names

    ## Filter on attributes
    if( $passes_filter && $attributes ) {
      foreach my $tag ( keys %$attributes ) {
        if( !$feature->has_tag( $tag ) ) {
          $passes_filter = 0;
          last;
        }
        $passes_filter = 0;
        foreach my $value ( $feature->get_tag_values( $tag ) ) {
          if( $value eq $attributes->{ $tag } ) {
            $passes_filter = 1;
            last;
          }
        } # End for each $value
        last unless( $passes_filter );
      } # End for each $tag
    } # End if we're filtering on attributes
    
    ## Filter on range
    if( $passes_filter && $ranges ) {
      foreach my $range ( @$ranges ) {
        if( $range->seq_id() &&
            $feature->seq_id() &&
            !( $feature->seq_id() eq $range->seq_id() )
          ) {
          ## If they both define seq_id(), and they don't match, then
          ## the ranges don't match.
          $passes_filter = 0;
        } elsif(
            ( $rangetype eq 'overlaps' ) &&
            !$range->overlaps( $feature, $strandmatch )
          ) {
          $passes_filter = 0;
        } elsif(
            ( $rangetype eq 'contains' ) &&
            !$range->contains( $feature, $strandmatch )
          ) {
          $passes_filter = 0;
        } elsif(
            ( $rangetype eq 'contained_in' ) &&
            !$feature->contains( $range, $strandmatch )
          ) {
          $passes_filter = 0;
        }
        last unless $passes_filter;
      } # End foreach $range in @$ranges
    } # End if we're filtering on range

    next unless( $passes_filter );

    if( defined $callback ) {
      unless( $callback->( $feature, $self ) ) {
        return 0; # return false when the callback fails.
      }
    } else {
      push( @features_to_return, $feature );
    }
  } # End for each $feature

  if( $iterator ) {
    return Bio::SeqFeature::SimpleIterator->new( @features_to_return );
  } elsif( defined $callback ) {
    return 1;
  } else {
    return @features_to_return;
  }
} # features()

=head2 feature_count

 Title   : feature_count
 Usage   : $collection->feature_count()
 Function: Return the number of SeqFeatureI objects that would be returned by a
           call to features() with no arguments.
 Returns : integer representing the number of SeqFeatureI objects
 Args    : None

=cut

sub feature_count {
  my $self = shift;

  my $count = 0;

  foreach my $seq_id ( keys %{ $self->{ '_seq_id_to_feature_table' } } ) {
    $count += scalar keys %{ $self->{ '_seq_id_to_feature_table' }->
                                    { $seq_id }->
                                    { '_identifiable_features' } };
    foreach my $start_position (
      keys %{ $self->{ '_seq_id_to_feature_table' }->
                     { $seq_id }->
                     { '_anonymous_features' } }
    ) {
      $count +=
        scalar @{ $self->{ '_seq_id_to_feature_table' }->
                         { $seq_id }->
                         { '_anonymous_features' }->
                         { $start_position } };
    }
  }

  return $count;
} # feature_count()

=head2 _insert_feature

 Title   : _insert_feature
  Usage   : $provider->_insert_feature( $feature );
 Function: Inserts the given feature into the store.
 Args    : L<Bio::SeqFeature::SegmentI> object
 Returns : The feature added, or undef iff the feature already existed.
 Status  : Protected

  Note that the returned feature might not be the same object as the
  given feature if, for example, the _create_feature(..) method was
  invoked on the argument.

=cut

# Enveloped to add $self as an observer of the added feature.
sub _insert_feature {
  my $self = shift;
  my $added_feature = $self->SUPER::_insert_feature( @_ );
  if( defined $added_feature ) {
    $added_feature->add_observer( $self );
  }
  return $added_feature;
} # _insert_feature(..)

=head2 _insert_or_update_feature

 Title   : _insert_or_update_feature
 Usage   : $provider->_insert_or_update_feature( $feature );
 Function: Inserts or updates the given feature in the store.
 Args    : L<Bio::SeqFeature::SegmentI> object
 Returns : The feature that was added or updated.
 Status  : Protected

  Note that the returned feature might not be the same object as the
  given feature if, for example, the _create_feature(..) method was
  invoked on the argument.

=cut

# Enveloped to ensure that $self is an observer of the feature.
sub _insert_or_update_feature {
  my $self = shift;
  my $inserted_or_updated_feature =
    $self->SUPER::_insert_or_update_feature( @_ );

  if( defined $inserted_or_updated_feature ) {
    $inserted_or_updated_feature->delete_observer( $self );
    $inserted_or_updated_feature->add_observer( $self );
  }
  return $inserted_or_updated_feature;
} # _insert_or_update_feature(..)

=head2 _remove_feature

 Title   : _remove_feature
 Usage   : $provider->_remove_feature( $feature );
 Function: Removes the given feature from the store.
 Args    : L<Bio::SeqFeature::SegmentI> object
 Returns : The removed feature, or false iff the feature was not
           previously in the store.
 Status  : Protected

=cut

# Enveloped to remove $self as an observer from the removed feature.
sub _remove_feature {
  my $self = shift;
  my $removed_feature = $self->SUPER::_remove_feature( @_ );
  if( defined $removed_feature ) {
    $removed_feature->delete_observer( $self );
  }
  return $removed_feature;
} # _remove_feature(..)

# For observing features.  Their notify_observers(..) method calls this.
sub update {
  my ( $self, $object, $action, %params ) = @_;
  # We need to know if the features stored by start position change
  # their start position, and if the features stored by unique_id
  # change their unique_id, and if features stored by start position
  # gain a unique_id.
  ## TODO: REMOVE
  #print STDERR "for action '$action', params are { ", join( ", ", %params ), " }.\n";
  return unless( $action );
  if( $action eq 'start' ) {
    return if( defined( $object->unique_id() ) );
    ## TODO: REMOVE
    #print STDERR "start changed: ",$params{ '-old' }," to ",$params{ '-new' },"\n";
    ## TODO: Dehackify.
    $self->remove_features( $object );
    $self->add_features( $object );
  } elsif( $action eq 'unique_id' ) {
    ## TODO: REMOVE
    #print STDERR "unique_id changed: ",$params{ '-old' }," to ",$params{ '-new' },"\n";
    ## TODO: Dehackify.
    $self->remove_features( $object );
    $self->add_features( $object );
  }
} # update(..)

1;

__END__
