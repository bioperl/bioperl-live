package Bio::SeqFeature::SimpleCollectionProvider;

# $Id $
## A simple implementation of the Bio::SeqFeature::CollectionProviderI
## interface.

=head1 NAME

Bio::SeqFeature::SimpleCollectionProvider - An in-memory
implementation of the CollectionProviderI interface.

=head1 SYNOPSIS

  my $provider = new Bio::SeqFeature::SimpleCollectionProvider();
  my $fg_color = $provider->get( 'fgcolor' );

=head1 DESCRIPTION

A SimpleCollectionProvider is an in-memory store of SeqFeatureIs that
implements the CollectionProviderI interface.  It supports updating,
adding, and removing SeqFeatures in its store, although explicit
updating is unnecessary as the SeqFeatureIs returned by the
get_collection method are the same as those stored herein (ie. there
is no external backing store).

Features can be filtered by the following attributes:

  1) their location, perhaps relative to a range (with a choice
     between overlapping, contained within, or completely containing a
     range)

  2) their type

  3) other attributes using tag/value semantics

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
use Bio::RelRange qw( &absSeqId );

use Bio::Root::Root;
use Bio::SeqFeature::CollectionProviderI;

$VERSION = '0.01';
@ISA = qw( Bio::Root::Root Bio::SeqFeature::CollectionProviderI );

=head2 new

 Title   : new
 Usage   : my $obj =
             new Bio::SeqFeature::SimpleCollectionProvider( @features );
 Function: Builds a new Bio::SeqFeature::SimpleCollectionProvider object 
 Returns : Bio::SeqFeature::SimpleCollectionProvider
 Args    : SeqFeatureI objects to store herein

=cut

sub new {
  my( $class, @args ) = @_;

  my $self = $class->SUPER::new( @args );
  $self->_initialize_simple_collection_provider();
  return $self;
} # new(..)

sub _initialize_simple_collection_provider {
  my $self = shift;
  my @args = @_;

  return if( $self->{ '_simple_collection_provider_initialized' } );

  $self->{ '_seq_id_to_feature_table' } = {};

  my ( $features );
  if( scalar( @args ) && ( $args[ 0 ] =~ /^-/ ) ) {
    ( $features ) =
      rearrange(
        [
         [ qw( FEATURES FEATURE ) ],
        ],
        @args
      );
  } else {
    # features
    $features = \@args;
  }
  ## Fix up features.
  if( $features ) {
    unless( ref $features eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $features = [ $features ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$features ) ||
       ( ( scalar( @$features ) == 1 ) && !( $features->[ 0 ] ) )
      ) {
      undef $features;
    }
  }

  if( $features ) {
    foreach my $feature ( @$features ) {
      unless( ref( $feature ) && $feature->isa( "Bio::SeqFeature::SegmentI" ) ) {
        $self->throw( "The given feature $feature is not a Bio::SeqFeature::SegmentI.  It is a ".ref( $feature ) );
      }
      unless( $self->_insert_feature( $feature ) ) {
        $self->throw( "duplicate feature: $feature" );
      }
    }
  }

  $self->{ '_simple_collection_provider_initialized' }++;
  return $self;
} # _initialize_simple_collection_provider(..)

=head2 new_from_collectionprovider

 Title   : new_from_collectionprovider
 Usage   : my $new_collection_provider =
             Bio::SeqFeature::SimpleCollectionProvider->
               new_from_collectionprovider(
                 $copy_from
               );
 Function: Create a new Bio::SeqFeature::SimpleCollectionProvider object
           by copying values from another SimpleCollectionProvider object.
 Returns : A new L<Bio::SeqFeature::SimpleCollectionProvider> object
 Args    : Another L<Bio::SeqFeature::SimpleCollectionProvider> object
 Status  : Protected

  This is a special copy constructor.  It forces the new collection
  provider into the L<Bio::SeqFeature::SimpleCollectionProvider>
  package, regardless of the package that it is called from.  This
  causes subclass-specfic information to be dropped.

  This also does not copy into the new provider the features held in
  the existing provider.  If you would like the new provider to hold the
  same features you must explicitly add them, like so:
    $new_collection_provider->insert_collection(
      $copy_from->get_collection()
    );

  As a special bonus you may also pass an existing hash and it will be the
  blessed an anointed object that is returned, like so:
    $new_collection_provider =
      Bio::SeqFeature::SimpleCollectionProvider->new_from_collectionprovider(
        $copy_from,
        $new_collection_provider
      );

=cut

sub new_from_collectionprovider {
  my $pack = shift; # ignored
  my $provider = shift || $pack;
  my $new_provider = shift || Bio::SeqFeature::SimpleCollectionProvider->new();
  @{ $new_provider }{ qw( _parent_collection_provider ) } =
    @{ $provider }{ qw( _parent_collection_provider ) };
  ## Note that we do not copy the features.
  return bless $new_provider, __PACKAGE__;
} # new_from_collection(..)

=head2 get_collection

 Title   : get_collection
 Usage   : my $collection = $collectionprovider->get_collection( @args );
 Returns : A L<Bio::SeqFeature::CollectionI> object
 Args    : see below
 Status  : Public

This routine will retrieve a L<Bio::SeqFeature::CollectionI> object
based on feature type, location or attributes.  The SeqFeatureI
objects in the returned CollectionI may or may not be newly
instantiated by this request.  If you make a modification to a feature
you must call update_collection with a collection that contains that
feature to ensure that the data provider is in sync with your change.
You may not, however, assume that modifications to the feature do not
auto-sync (they might!).  This simple implementation will auto-sync,
although subclasses may not.

If a range is specified using the -range argument then this range will
 be used to narrow the results, according to the specified -rangetype
 and -strandtype arguments.

-rangetype is one of:
   "overlaps"      return all features that overlap the range (default)
   "contains"      return features completely contained within the range
   "contained_in"  return features that completely contain the range

Note that if the implementing class implements L<Bio::LocationI> then the
baselocation will default to that range.  If a baselocation is given
or defaulted and a range is specified as an argument, then the
coordinates of the given range will be interpreted relative to the
implementing class\'s range.  If the implementing class does not
implement LocationI and no range is given, then -rangetype may be
ignored.

-strandmatch is one of:
   "strong"        ranges must have the same strand
   "weak"          ranges must have the same strand or no strand (default)
   "ignore"        ignore strand information

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of feature types (as if they
were given as -types => \@args).  In the named parameter form, the
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

  -location      A L<Bio::LocationI> object defining the range to search and
                 the rangetype.  Use -range (and -baselocation,
                 perhaps; see below) as an alternative to -location.
                 See also -strandmatch.  There may be a default value
                 for -location.

  -baselocation  A L<Bio::LocationI> object defining the location to which
                 the -range argument is relative.  There may be a
                 default -baselocation.  If this CollectionProviderI is also a
                 L<Bio::LocationI>, then the default -baselocation should be
                 itself.

  -range         A L<Bio::RangeI> object defining the range to search.
                 See also -strandmatch and -rangetype.  Use instead of
                 -location, when -baselocation is specified or
                 provided by default (see above).

  -rangetype     One of "overlaps", "contains", or "contained_in".

  -strandmatch   One of "strong", "weak", or "ignore".

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

## This is a hacky implementation that delegates all the work to a
## CollectionI implementation (if this *is* one, then that's
## delegating to our own features() method).  This works because the
## argument list of get_collection is a subset of the argument list of
## features().
## This implementation uses _create_collection() for all collection
## construction.  See it for info about how to override it.
sub get_collection {
  my $self = shift;

  ## TODO: REMOVE
  #warn "SimpleCollectionProvider::get_collection( ( ".join( ', ', @_ )." )";# if Bio::Graphics::Browser::DEBUG;

  if( $self->isa( 'Bio::SeqFeature::CollectionI' ) ) {
    ## Use our own features() method to do the hard work..
    return $self->_create_collection( \@_, 'lookup' );
  } else {
    # We're not a CollectionI.  Let's hijack one for our own nefarious needs.
    my $hijacked_collection = $self->_create_collection( \@_ );
    # Add everything to it.
    foreach my $seq_id ( keys %{ $self->{ '_seq_id_to_feature_table' } } ) {
      $hijacked_collection->add_features(
        values %{ $self->{ '_seq_id_to_feature_table' }->
                         { $seq_id }->
                         { '_identifiable_features' }
                }
      );
      foreach my $start ( keys %{ $self->{ '_seq_id_to_feature_table' }->
                                         { $seq_id }->
                                         { '_anonymous_features' }
                                } ) {
        $hijacked_collection->add_features(
          @{ $self->{ '_seq_id_to_feature_table' }->
                    { $seq_id }->
                    { '_anonymous_features' }->
                    { $start }
           }
        );
      } # End foreach $start
    } # End foreach $seq_id

    # Now delegate to it.
    return $hijacked_collection->get_collection( @_ );
    ## A slightly more clever hack would keep the hijacked collection
    ## as its backing store so that future calls to get_collection are
    ## faster...
  }
} # get_collection(..)

=head2 insert_or_update_collection

 Title   : insert_or_update_collection
 Usage   : $collectionprovider->insert_or_update( $collection );
 Function: Attempts to update all the features of a collection.  If
           a feature doesn\'t exist it inserts it automatically.
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

sub insert_or_update_collection {
  my $self = shift;
  my $collection = shift;

  return unless defined( $collection );

  my $iterator = $collection->get_feature_stream();
  my $feature;
  while( $iterator->has_more_features() ) {
    $feature = $iterator->next_feature();
    $self->_insert_or_update_feature( $feature );
  }
} # insert_or_update_collection(..)

=head2 insert_collection

 Title   : insert_collection
 Usage   : $collectionprovider->insert_collection($collection);
 Function: Insert all the features of a collection.  If any features
           already exist throw an exception. 
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

sub insert_collection {
  my $self = shift;
  my $collection = shift;

  return unless defined( $collection );

  my $iterator = $collection->get_feature_stream();
  my $feature;
  while( $iterator->has_more_features() ) {
    $feature = $iterator->next_feature();
    unless( $self->_insert_feature( $feature ) ) {
      $self->throw( "duplicate feature: $feature" );
    }
  }
} # insert_collection(..)

=head2 update_collection

 Title   : update_collection
 Usage   : $collectionprovider->update_collection($collection);
 Function: Updates all the features of a collection.  If any do not
           already exist throw an exception.
 Returns : Return the updated collection upon success or undef
           upon failure.
 Args    : L<Bio::SeqFeature::CollectionI> object

  If you make a modification to a feature you must call
  update_collection with a collection that contains that feature to
  ensure that the data provider is in sync with your change.  You may
  not, however, assume that modifications to the feature do not
  auto-sync (they might!).

=cut

sub update_collection {
  my $self = shift;
  my $collection = shift;

  return unless defined( $collection );

  my $iterator = $collection->get_feature_stream();
  my $feature;
  while( $iterator->has_more_features() ) {
    $feature = $iterator->next_feature();
    unless( $self->_update_feature( $feature ) ) {
      $self->throw( "nonexistent feature: $feature" );
    }
  }
} # update_collection(..)

=head2 remove_collection

 Title   : remove_collection
 Usage   : $provider->remove_collection($collection);
 Function: Removes all the features in a collection.  If any features 
           do not exists throw an exception.
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

sub remove_collection {
  my $self = shift;
  my $collection = shift;

  return unless defined( $collection );

  my $iterator = $collection->get_feature_stream();
  my $feature;
  while( $iterator->has_more_features() ) {
    $feature = $iterator->next_feature();
    unless( $self->_remove_feature( $feature ) ) {
      $self->throw( "nonexistent feature: $feature" );
    }
  }
} # remove_collection(..)

=head2 types

 Title   : types
 Usage   : my @types = $collectionprovider->types();
           OR
           my %types_and_counts = $collectionprovider->types( -count => 1 );
 Function: Enumerate the feature types provided by this provider, and possibly
           count the features in each type.
 Returns : a list of L<Bio::SeqFeature::TypeI> objects
           OR
           a hash mapping type id strings to integer counts
 Args    : see below

This routine returns a list of feature types known to the provider.
If the -count argument is given, it returns a hash of known types
mapped to their occurrence counts in this provider.  Note that the
returned list (or the keys of the returned hash) may include types for
which the count is 0.  Also note that the hierarchy of TypeI objects
is ignored, so if there are 4 features of type 'foo' which is a child
of type 'bar', and only 1 feature (explicitly) of type 'bar', then the
count for 'bar' will be 1, not 5.

Arguments are -option=E<gt>value pairs as follows:

  -count aka -enumerate  if true, count the features

The returned value will be a list of L<Bio::SeqFeature::TypeI> objects
or a hash with the string values of these objects as keys.

The SimpleCollectionProvider implementation recalculates the types
with each call to types(), so the method takes time proportional to
the number of features provided by this provider.  Subclasses may not
behave this way; they should document the difference if there is one.

=cut

sub types {
  my $self = shift;
  my $count;
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $count ) =
      rearrange( [ qw( COUNT COUNTS ENUM ENUMERATE ) ], @_ );
  }
  my %types;
  foreach my $seq_id ( keys %{ $self->{ '_seq_id_to_feature_table' } } ) {
    foreach my $feature (
      ( values %{ $self->{ '_seq_id_to_feature_table' }->
                         { $seq_id }->
                         { '_identifiable_features' }
                } ),
      ( map @{ $_ },
        ( values %{ $self->{ '_seq_id_to_feature_table' }->
                           { $seq_id }->
                           { '_anonymous_features' }
                  } ) )
    ) {
      $types{ $feature->type() }++;
    }
  }
  if( $count ) {
    return %types;
  } else {
    return keys %types;
  }
} # types(..)

sub seq_ids {
  my $self = shift;
  my $count;
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $count ) =
      rearrange( [ qw( COUNT COUNTS ENUM ENUMERATE ) ], @_ );
  }
  # Note here that the seq_ids that are used to store features are
  # absolute seq_ids, so we can just use our _seq_id_to_feature_table hash.
  my %seq_ids;
  foreach my $seq_id ( keys %{ $self->{ '_seq_id_to_feature_table' } } ) {
    ## TODO: There's probably a more efficient way to do this.
    foreach my $feature (
      ( values %{ $self->{ '_seq_id_to_feature_table' }->
                         { $seq_id }->
                         { '_identifiable_features' }
                } ),
      ( map @{ $_ },
        ( values %{ $self->{ '_seq_id_to_feature_table' }->
                           { $seq_id }->
                           { '_anonymous_features' }
                  } ) )
    ) {
      $seq_ids{ $seq_id }++;
    }
  }
  if( $count ) {
    return %seq_ids;
  } else {
    return keys %seq_ids;
  }
} # seq_ids(..)

=head2 parent_collection_provider

 Title   : parent_collection_provider
 Usage   : my $parent = $provider->parent_collection_provider( [$new_parent] );
 Function: Get/set the CollectionProviderI that is the parent of this provider.
 Returns : The current (or former, if used as a set method)
           L<Bio::SeqFeature::CollectionProviderI> parent or undef if
           there is none
 Args    : [optional] a new parent

  CollectionProviderIs may be views onto other CollectionProviderIs.
  A common example is the CollectionI returned by the get_collection()
  method.  It is a CollectionProviderI as well (CollectionI ISA
  CollectionProviderI), but it (initially) provides only features
  found in the original CollectionProviderI.  The original is then
  called its parent, and is returned by calling this method.  Note the
  following facts:

    1) Not all CollectionProviderIs have parents.

    2) A CollectionProviderI may store its features independently from
       its parent or it may delegate to its parent; the behavior is
       unspecified by the interface.

    3) A CollectionProviderI may have features that its parent does
       not have; this may happen eg. when a feature was added to the
       CollectionProviderI but not to its parent.

=cut

sub parent_collection_provider {
  my $self = shift;
  my $new_value = shift;
  my $old_value = $self->{ '_parent_collection_provider' };
  if( defined( $new_value ) ) {
    unless( ref( $new_value ) && $new_value->isa( 'Bio::SeqFeature::CollectionProviderI' ) ) {
      $self->throw( "The given value is illegal because it is not a Bio::SeqFeature::CollectionProviderI object.  It is $new_value, a ".ref( $new_value )."." );
    }
    $self->{ '_parent_collection_provider' } = $new_value;
  }
  return $old_value;
} # parent_collection_provider(..)

=head2 _create_collection

 Title   : _create_collection
 Usage   : my $collection = $provider->_create_collection(
                           \@args_to_get_collection,
                           @features
                         );
           or
           my $collection = $provider->_create_collection(
                           \@args_to_get_collection,
                           'lookup'
                         );
 Function: Factory method for instantiating a collection.
 Args    : a ref to the args used by the get_collection method, and some
           L<Bio::SeqFeatureI> objects to add to the new collection.
 Returns : a new L<Bio::SeqFeature::CollectionI> object
 Status  : Protected

   NOTE THIS CONSTRAINT: Because of our hacky implementation of
get_collection(), we require that the L<Bio::SeqFeature::CollectionI>
that is returned by this method creates collections of its own type
when its get_collection() method is called.  This is the likely
behavior anyway, but it should be noted, just in case.

=cut

# This implementation just ignores the args and makes a new
# SimpleCollection with the given features.
sub _create_collection {
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
  #warn "SimpleCollectionProvider::_create_segment( { ".join( ', ', ( my @foo = %$args ) )." }, ( ".join( ', ', @features )." )" if Bio::Graphics::Browser::DEBUG;

  return Bio::SeqFeature::SimpleCollection->new( @features );
} # _create_collection

=head2 _create_feature

 Title   : create_feature
 Usage   : my $new_seq_feature = $provider->_create_feature( $feature_data );
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

## This one makes Generic features when given [ $start, $end ] pairs.
sub _create_feature {
  my $self = shift;
  my $feature_data = shift;
  if( ref( $feature_data ) eq 'ARRAY' ) {
    # This is for feature data like [ [ start0, end0 ], [ start1, end1 ] ].
    my ( $start, $end ) = @{ $feature_data };
    
    # The following line should be unnecessary.
    #next unless defined $start && defined $end;
    
    # Strandedness defaults to 1, but start > end forces -1.
    my $strand = 1;
    if( $start > $end ) {
      ( $start, $end ) = ( $end, $start );
      $strand = -1;
    }
    return $self->new( -start  => $start,
                       -end    => $end,
                       -strand => $strand,
                       -parent => $self ) || $feature_data;
  }
  # If we're at this point then we've been unable to make a new feature.
  return $feature_data;
} # _create_feature(..)

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

sub _insert_feature {
  my $self = shift;
  my $feature = shift;

  unless( defined( $feature ) ) {
    $self->throw( "\$simple_collection_provider->_insert_feature( undef ): \$feature is undef!" );
  }
  unless( ref( $feature ) &&
          ( ref( $feature ) ne 'ARRAY' ) &&
          $feature->isa( "Bio::SeqFeature::SegmentI" ) ) {
    # Try to make a feature out of it.
    $feature = $self->_create_feature( $feature );
  }
  unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeature::SegmentI' ) ) {
    $self->throw( "\$simple_collection_provider->_insert_feature( $feature ): \$feature is not a Bio::SeqFeature::SegmentI!" );
  }

  my $seq_id = absSeqId( $feature );
  unless( defined $seq_id ) {
    $seq_id = 'undef';
  }
  if( defined( $feature->unique_id() ) ) {
    if( $self->{ '_seq_id_to_feature_table' }->
               { $seq_id }->
               { '_identifiable_features' }->
               { $feature->unique_id() }
      ) {
      return undef;
    } else {
      $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_identifiable_features' }->
             { $feature->unique_id() } =
               $feature;
      return $feature;
    }
  } else { # it does not have a unique id.  Store it by start location.
    my $features_that_start_where_this_one_does =
      $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_anonymous_features' }->
             { $feature->start() };
    if( $features_that_start_where_this_one_does &&
        scalar( @$features_that_start_where_this_one_does ) ) {
      foreach my $other_feature ( @$features_that_start_where_this_one_does ) {
        if( $other_feature == $feature ) {
          ## TODO: REMOVE. Testing.
          #warn "SimpleCollectionProvider::_insert_feature(..): $feature == $other_feature.  \$feature isa ".ref( $feature ).", \$other_feature isa ".ref( $other_feature ).".  \$feature's display_name is ".$feature->display_name().", \$other_feature's display_name is ".$other_feature->display_name()."\n";
          return undef;
        }
      }
      push( @$features_that_start_where_this_one_does, $feature );
      return $feature;
    } else {
      $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_anonymous_features' }->
             { $feature->start() } = [ $feature ];
      return $feature;
    }
  }
} # _insert_feature(..)

=head2 _update_feature

 Title   : _update_feature
 Usage   : $provider->_update_feature( $feature );
 Function: Updates the given feature in the store.
 Args    : L<Bio::SeqFeature::SegmentI> object
 Returns : The updated feature, or false iff the feature is not in the
           store (it won\'t be added!)
 Status  : Protected

=cut

sub _update_feature {
  my $self = shift;
  my $feature = shift;

  unless( defined( $feature ) ) {
    $self->throw( "\$simple_collection_provider->_update_feature( undef ): \$feature is undef!" );
  }
  unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeature::SegmentI' ) ) {
    $self->throw( "\$simple_collection_provider->_update_feature( $feature ): \$feature is not a Bio::SeqFeature::SegmentI!" );
  }

  my $seq_id = absSeqId( $feature );
  unless( defined $seq_id ) {
    $seq_id = 'undef';
  }
  if( defined( $feature->unique_id() ) ) {
    if( $self->{ '_seq_id_to_feature_table' }->
               { $seq_id }->
               { '_identifiable_features' }->
               { $feature->unique_id() }
      ) {
      $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_identifiable_features' }->
             { $feature->unique_id() } =
        $feature;
      return $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_identifiable_features' }->
             { $feature->unique_id() };
    } else {
      return undef;
    }
  } else { # it does not have a unique id.  It is stored by start location.
    my $features_that_start_where_this_one_does =
      $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_anonymous_features' }->
             { $feature->start() };
    if( $features_that_start_where_this_one_does &&
        scalar( @$features_that_start_where_this_one_does ) ) {
      for( my $i = 0;
           $i < scalar( @$features_that_start_where_this_one_does );
           $i++
         ) {
        if(
           $features_that_start_where_this_one_does->[ $i ] == $feature
          ) {
          $features_that_start_where_this_one_does->[ $i ] = $feature;
          return $features_that_start_where_this_one_does->[ $i ];
        }
      }
      return undef;
    } else {
      return undef;
    }
  }
} # _update_feature(..)

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

sub _insert_or_update_feature {
  my $self = shift;
  my $feature = shift;

  unless( defined( $feature ) ) {
    $self->throw( "\$simple_collection_provider->_insert_or_update_feature( undef ): \$feature is undef!" );
  }
  unless( ref( $feature ) &&
          ( ref( $feature ) ne 'ARRAY' ) &&
          $feature->isa( "Bio::SeqFeature::SegmentI" ) ) {
    # Try to make a feature out of it.
    $feature = $self->_create_feature( $feature );
  }
  unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeature::SegmentI' ) ) {
    $self->throw( "\$simple_collection_provider->_insert_or_update_feature( $feature ): \$feature is not a Bio::SeqFeature::SegmentI!" );
  }

  my $seq_id = absSeqId( $feature );
  unless( defined $seq_id ) {
    $seq_id = 'undef';
  }
  if( defined( $feature->unique_id() ) ) {
    $self->{ '_seq_id_to_feature_table' }->
           { $seq_id }->
           { '_identifiable_features' }->
           { $feature->unique_id() } =
      $feature;
    return $self->{ '_seq_id_to_feature_table' }->
           { $seq_id }->
           { '_identifiable_features' }->
           { $feature->unique_id() };
  } else { # it does not have a unique id.  It is stored by start location.
    my $features_that_start_where_this_one_does =
      $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_anonymous_features' }->
             { $feature->start() };
    if( $features_that_start_where_this_one_does &&
        scalar( @$features_that_start_where_this_one_does ) ) {
      for( my $i = 0;
           $i < scalar( @$features_that_start_where_this_one_does );
           $i++
         ) {
        if(
           $features_that_start_where_this_one_does->[ $i ] == $feature
          ) {
          $features_that_start_where_this_one_does->[ $i ] = $feature;
          return $features_that_start_where_this_one_does->[ $i ];
        }
      }
      push( @$features_that_start_where_this_one_does, $feature );
      return $feature;
    } else {
      $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_anonymous_features' }->
             { $feature->start() } = [ $feature ];
      return $feature;
    }
  }
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

sub _remove_feature {
  my $self = shift;
  my $feature = shift;

  unless( defined( $feature ) ) {
    $self->throw( "\$simple_collection_provider->_remove_feature( undef ): \$feature is undef!" );
  }
  unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeature::SegmentI' ) ) {
    $self->throw( "\$simple_collection_provider->_remove_feature( $feature ): \$feature is not a Bio::SeqFeature::SegmentI!" );
  }

  my $seq_id = absSeqId( $feature );
  unless( defined $seq_id ) {
    $seq_id = 'undef';
  }
  if( defined( $feature->unique_id() ) ) {
    if( $self->{ '_seq_id_to_feature_table' }->
               { $seq_id }->
               { '_identifiable_features' }->
               { $feature->unique_id() }
      ) {
      delete $self->{ '_seq_id_to_feature_table' }->
                    { $seq_id }->
                    { '_identifiable_features' }->
                    { $feature->unique_id() };
      return $self->{ '_seq_id_to_feature_table' }->
                    { $seq_id }->
                    { '_identifiable_features' }->
                    { $feature->unique_id() };
    } else {
      return undef;
    }
  } else { # it does not have a unique id.  It is stored by start location.
    my $features_that_start_where_this_one_does =
      $self->{ '_seq_id_to_feature_table' }->
             { $seq_id }->
             { '_anonymous_features' }->
             { $feature->start() };
    if( $features_that_start_where_this_one_does &&
        scalar( @$features_that_start_where_this_one_does ) ) {
      for( my $i = 0;
           $i < scalar( @$features_that_start_where_this_one_does );
           $i++
         ) {
        if(
           $features_that_start_where_this_one_does->[ $i ] == $feature
          ) {
          if( scalar( @$features_that_start_where_this_one_does ) == 1 ) {
            delete $self->{ '_seq_id_to_feature_table' }->
                          { $seq_id }->
                          { '_anonymous_features' }->
                          { $feature->start() };
          } else {
            splice( @$features_that_start_where_this_one_does, $i, 1 );
          }
          return $features_that_start_where_this_one_does->[ $i ];
        }
      }
      return undef;
    } else {
      return undef;
    }
  }
} # _remove_feature(..)

=head2 toString

 Title   : toString
 Usage   : $str_val = $collection->toString()
 Function: returns "A SimpleCollectionProvider.";
 Returns : a String
 Args    : None
 Status  : Public

  This method is a hack.

=cut

sub toString {
  my $self = shift;

  ## TODO: Dehackify
  if( $self->isa( "Bio::SeqFeatureI" ) ) {
    return Bio::SeqFeatureI::toString( $self );
  }
  return "A SimpleCollectionProvider.";
} # toString()

## method for overload for comparing two SimpleCollectionProvider objects.  Uses toString().
sub _cmp {
  my $self = shift;
  my ( $b, $reversed ) = @_;
  my $a = $self->toString();
  ( $a, $b ) = ( $b, $a ) if $reversed;
  return ( $a cmp $b );
}

1;

__END__
