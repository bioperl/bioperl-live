package Bio::SeqFeature::SimpleCollectionProvider;

# $Id $

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
  $self->{ '_identifiable_features' } = {};
  $self->{ '_anonymous_features' } = {};

  if( scalar( @args ) ) {
    foreach my $feature ( @args ) {
      next unless( ref( $feature ) && $feature->isa( "Bio::SeqFeatureI" ) );
      unless( $self->_insert_feature( $feature ) ) {
        $self->throw( "duplicate feature: $feature" );
      }
    }
  }
  return $self;
} # new(..)

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
sub get_collection {
  my $self = shift;
  if( $self->isa( 'Bio::SeqFeature::CollectionI' ) ) {
    ## Use our own features() method to do the hard work..
    return new Bio::SeqFeature::SimpleCollection( $self->features( @_ ) );
  } else {
    # We're not a CollectionI.  Let's hijack one for our own nefarious needs.
    my $hijacked_collection = new Bio::SeqFeature::SimpleCollection();
    # Add everything to it.
    $hijacked_collection->add_features(
      values %{ $self->{ '_identifiable_features' } }
    );
    foreach my $start ( keys %{ $self->{ '_anonymous_features' } } ) {
      $hijacked_collection->add_features(
        @{ $self->{ '_anonymous_features' }{ $start } }
      );
    }

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
  foreach my $feature ( ( values %{ $self->{ '_identifiable_features' } } ),
                        ( map +( @{ $_ }, 1 ),
                          ( values %{ $self->{ '_anonymous_features' } } ) ) ) {
    $types{ $feature->type() }++;
  }
  if( $count ) {
    return %types;
  } else {
    return keys %types;
  }
} # types(..)

=head2 _insert_feature

 Title   : _insert_feature
  Usage   : $provider->_insert_feature( $feature );
 Function: Inserts the given feature into the store.
 Args    : L<Bio::SeqFeatureI> object
 Returns : False iff the feature already existed.
 Status  : Protected

=cut

sub _insert_feature {
  my $self = shift;
  my $feature = shift;

  unless( defined( $feature ) ) {
    $self->throw( "\$simple_collection_provider->_insert_feature( undef ): \$feature is undef!" );
  }
  unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeatureI' ) ) {
    $self->throw( "\$simple_collection_provider->_insert_feature( $feature ): \$feature is not a Bio::SeqFeatureI!" );
  }

  if( defined( $feature->unique_id() ) ) {
    if( $self->{ '_identifiable_features' }{ $feature->unique_id() } ) {
      return 0;
    } else {
      $self->{ '_identifiable_features' }{ $feature->unique_id() } =
        $feature;
      return 1;
    }
  } else { # it does not have a unique id.  Store it by start location.
    my $features_that_start_where_this_one_does =
      $self->{ '_anonymous_features' }{ $feature->start() };
    if( $features_that_start_where_this_one_does &&
        scalar( @$features_that_start_where_this_one_does ) ) {
      foreach my $other_feature ( @$features_that_start_where_this_one_does ) {
        if( ( $other_feature == $feature )
            ||
            $other_feature->equals( $feature )
          ) {
          return 0;
        }
      }
      push( @$features_that_start_where_this_one_does, $feature );
      return 1;
    } else {
      $self->{ '_anonymous_features' }{ $feature->start() } = [ $feature ];
      return 1;
    }
  }
} # _insert_feature(..)

=head2 _update_feature

 Title   : _update_feature
 Usage   : $provider->_update_feature( $feature );
 Function: Updates the given feature in the store.
 Args    : L<Bio::SeqFeatureI> object
 Returns : False iff the feature is not in the store (it won\'t be added!)
 Status  : Protected

=cut

sub _update_feature {
  my $self = shift;
  my $feature = shift;

  unless( defined( $feature ) ) {
    $self->throw( "\$simple_collection_provider->_update_feature( undef ): \$feature is undef!" );
  }
  unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeatureI' ) ) {
    $self->throw( "\$simple_collection_provider->_update_feature( $feature ): \$feature is not a Bio::SeqFeatureI!" );
  }

  if( defined( $feature->unique_id() ) ) {
    if( $self->{ '_identifiable_features' }{ $feature->unique_id() } ) {
      $self->{ '_identifiable_features' }{ $feature->unique_id() } =
        $feature;
      return 1;
    } else {
      return 0;
    }
  } else { # it does not have a unique id.  It is stored by start location.
    my $features_that_start_where_this_one_does =
      $self->{ '_anonymous_features' }{ $feature->start() };
    if( $features_that_start_where_this_one_does &&
        scalar( @$features_that_start_where_this_one_does ) ) {
      my $other_feature;
      for( my $i = 0;
           $i < scalar( @$features_that_start_where_this_one_does );
           $i++
         ) {
        if(
           ( $features_that_start_where_this_one_does->[ $i ] == $feature )
           ||
           $features_that_start_where_this_one_does->[ $i ]->equals( $feature )
          ) {
          $features_that_start_where_this_one_does = $feature;
          return 1;
        }
      }
      return 0;
    } else {
      return 0;
    }
  }
} # _update_feature(..)

=head2 _insert_or_update_feature

 Title   : _insert_or_update_feature
 Usage   : $provider->_insert_or_update_feature( $feature );
 Function: Inserts or updates the given feature in the store.
 Args    : L<Bio::SeqFeatureI> object
 Returns : True
 Status  : Protected

=cut

sub _insert_or_update_feature {
  my $self = shift;
  my $feature = shift;

  unless( defined( $feature ) ) {
    $self->throw( "\$simple_collection_provider->_insert_or_update_feature( undef ): \$feature is undef!" );
  }
  unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeatureI' ) ) {
    $self->throw( "\$simple_collection_provider->_insert_or_update_feature( $feature ): \$feature is not a Bio::SeqFeatureI!" );
  }

  if( defined( $feature->unique_id() ) ) {
    $self->{ '_identifiable_features' }{ $feature->unique_id() } =
      $feature;
    return 1;
  } else { # it does not have a unique id.  It is stored by start location.
    my $features_that_start_where_this_one_does =
      $self->{ '_anonymous_features' }{ $feature->start() };
    if( $features_that_start_where_this_one_does &&
        scalar( @$features_that_start_where_this_one_does ) ) {
      my $other_feature;
      for( my $i = 0;
           $i < scalar( @$features_that_start_where_this_one_does );
           $i++
         ) {
        if(
           ( $features_that_start_where_this_one_does->[ $i ] == $feature )
           ||
           $features_that_start_where_this_one_does->[ $i ]->equals( $feature )
          ) {
          $features_that_start_where_this_one_does = $feature;
          return 1;
        }
      }
      push( @$features_that_start_where_this_one_does, $feature );
      return 1;
    } else {
      $self->{ '_anonymous_features' }{ $feature->start() } = [ $feature ];
      return 1;
    }
  }
} # _insert_or_update_feature(..)

=head2 _remove_feature

 Title   : _remove_feature
  Usage   : $provider->_remove_feature( $feature );
 Function: Removes the given feature from the store.
 Args    : L<Bio::SeqFeatureI> object
 Returns : False iff the feature was not previously in the store.
 Status  : Protected

=cut

sub _remove_feature {
  my $self = shift;
  my $feature = shift;

  unless( defined( $feature ) ) {
    $self->throw( "\$simple_collection_provider->_remove_feature( undef ): \$feature is undef!" );
  }
  unless( ref( $feature ) && $feature->isa( 'Bio::SeqFeatureI' ) ) {
    $self->throw( "\$simple_collection_provider->_remove_feature( $feature ): \$feature is not a Bio::SeqFeatureI!" );
  }

  if( defined( $feature->unique_id() ) ) {
    if( $self->{ '_identifiable_features' }{ $feature->unique_id() } ) {
      delete $self->{ '_identifiable_features' }{ $feature->unique_id() };
      return 1;
    } else {
      return 0;
    }
  } else { # it does not have a unique id.  It is stored by start location.
    my $features_that_start_where_this_one_does =
      $self->{ '_anonymous_features' }{ $feature->start() };
    if( $features_that_start_where_this_one_does &&
        scalar( @$features_that_start_where_this_one_does ) ) {
      my $other_feature;
      for( my $i = 0;
           $i < scalar( @$features_that_start_where_this_one_does );
           $i++
         ) {
        if(
           ( $features_that_start_where_this_one_does->[ $i ] == $feature )
           ||
           $features_that_start_where_this_one_does->[ $i ]->equals( $feature )
          ) {
          splice( @$features_that_start_where_this_one_does, $i, 1 );
          return 1;
        }
      }
      return 0;
    } else {
      return 0;
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
