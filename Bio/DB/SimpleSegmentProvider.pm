package Bio::DB::SimpleSegmentProvider;

# $Id$
# A simple implementation of Bio::DB::SegmentProviderI with an
# in-memory backing store, for testing purposes.

=head1 NAME

Bio::DB::SimpleSegmentProvider -- A simple provider of collections of sequence
segments, as if they were from a database or other non-trivial backing store.

=head1 SYNOPSIS

 use Bio::DB::SimpleSegmentProvider;

 use Bio::SeqFeature::Generic;
 use Bio::SeqFeature::SimpleCollection;

 my $data_provider =
   Bio::DB::SimpleSegmentProvider->new();

 # Add some features

 $data_provider->insert_collection(
   new Bio::SeqFeature::SimpleCollection(
     new Bio::SeqFeature::Generic(
       -id => 'foo',
       -start => 10,
       -end => 100,
       -strand => -1,
       -primary => 'repeat',
       -source_tag => 'repeatmasker',
       -score  => 1000
     ),
     new Bio::SeqFeature::Generic(
       -id => 'bar',
       -start => 100,
       -end => 200,
       -strand => -1
     )
   );
 );

 # Add another feature
 my $baz =
   new Bio::SeqFeature::Generic(
     -id => 'baz',
     -start => 1,
     -end => 200
   );
 $data_provider->insert_collection(
   new Bio::SeqFeature::SimpleCollection( $baz )
 );

 # Update one that we'd previously inserted.

 $baz->strand( -1 );
 $data_provider->update_collection(
   new Bio::SeqFeature::SimpleCollection( $baz );
 );

=head1 DESCRIPTION

The Bio::DB::SegmentProviderI interface provides access to
Bio::SeqFeature::CollectionIs and Bio::PrimarySeqIs stored in a
database or other (generally external) backing store.  It is a
Bio::DB::SimpleFeatureProvider and a Bio::DB::SimpleSequenceProvider.

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

use Bio::DB::SimpleSequenceProvider;
use Bio::DB::SimpleFeatureProvider;
use Bio::DB::SegmentProviderI;
@ISA = qw( Bio::DB::SimpleSequenceProvider
           Bio::DB::SimpleFeatureProvider
           Bio::DB::SegmentProviderI
         );

use vars '$VERSION';
$VERSION = '1.00';

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'

sub new {
  my( $class, @args ) = @_;

  my $self = $class->SUPER::new( @args );
  $self->_initialize_simple_segment_provider( @args );
  return $self;
} # new(..)

sub _initialize_simple_segment_provider {
  my $self = shift;
  my @args = @_;

  return if( $self->{ '_simple_segment_provider_initialized' } );

  ## NOTE: the sequence_provider initialization must happen first,
  ## because when we add features we try to add the sequences that
  ## they are on also, and to do this we have to have already
  ## initialized our sequence provision abilities.
  $self->_initialize_simple_sequence_provider( @args );
  $self->_initialize_simple_feature_provider( @args );

  $self->{ '_simple_segment_provider_initialized' }++;
  return $self;
} # _initialize_simple_segment_provider(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 get_collection

 Title    : get_collection
 Usage    : my $segment = $segmentprovider->get_collection( %args );
            OR
            my $segment = $segmentprovider->get_collection( @types );
            OR
            my @segments = $segmentprovider->get_collection( %args );
            OR
            my @segments = $segmentprovider->get_collection( @types );
 Returns  : A L<Bio::SeqFeature::SegmentI> object or a list thereof.
 Args     : see below
 Status   : Public
 Exception: "Those features do not share a common sequence" if this
            method is called in scalar context and the features that
            would otherwise be included in the resulting segment do
            not all fall on the same sequence.

NOTE: This method is (almost) identical to the get_collection() method
from L<Bio::DB::FeatureProviderI> that it overrides.  The entire
documentation follows, but first a brief summary of the changes:
  * This method returns L<Bio::SeqFeature::SegmentI> objects instead
    of mere CollectionI objects.  SegmentI objects are CollectionI
    objects, so this is an additional constraint on the interface.
    The returned SegmentI objects will have as their range the range
    searched, if any, or the smallest range that encloses the returned
    features.
  * This method will return a list of objects if called in list
    context; one L<Bio::SeqFeature::SegmentI> object per root sequence
    of the requested features.  Each returned SegmentI will have as
    its seq_id the common sequences' unique_id() or primary_id().
  * This method will throw an exception if called in scalar context
    and the features that would be included in the resulting SegmentI
    do not all share a common sequence.

This routine will retrieve one or more L<Bio::SeqFeature::SegmentI>
objects based on feature type, location or attributes.  The
SeqFeatureI objects in the returned SegmentIs may or may not be newly
instantiated by this request.  They will have as their range the range
searched, if any, or the smallest range that encloses the returned
features.  They will have as their seq_id() the unique_id() or
primary_id() of the returned features' common sequence.  If this
method is called in list context then one SegmentI object will be
returned per root sequence.  If this method is called in scalar
context and the returned features do not share a common sequence then
an exception will be thrown.

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
                 the -range argument is relative.  There may be a
                 default -baserange.  If this SegmentProviderI is also a
                 L<Bio::RangeI>, then the default -baserange should be
                 itself.

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

## This is a hacky implementation that delegates all the work to one or more
## SegmentI implementations (if this *is* one, then that's
## delegating to our own features() method).  This works because the
## argument list of get_collection is a subset of the argument list of
## features().
## This implementation uses _create_collection() for all
## Bio::SeqFeature::SegmentI construction.  See it for info about how
## to override it.
sub get_collection {
  my $self = shift;

  if( $self->isa( 'Bio::SeqFeature::SegmentI' ) ) {
    ## Use our own features() method to do the hard work.
    ## This is okay because if we are a SegmentI then all of our
    ## features are necessarily on the same sequence...
    return $self->_create_collection( \@_, 'lookup' );
  } else {
    # We're not a SegmentI.  We'll have to hijack one per seq_id for
    # our own nefarious needs.
    my $hijacked_segment;
    my @segments;
    my $sequence;
    # Add everything to it.
    foreach my $seq_id ( keys %{ $self->{ '_seq_id_to_feature_table' } } ) {
      # First see if we have that sequence in our stores..
      ( $sequence ) = $self->sequences( $seq_id );
      if( $sequence ) {
        ## TODO: REMOVE
        #warn "using stored sequence $sequence.\n" if DEBUG;
       } else {
        ## TODO: REMOVE
        #warn "we know nothing of the sequence with id $seq_id.\n" if DEBUG;
       }
      my %args =
        ( $sequence ? ( '-seq_id' => $sequence ) : ( '-seq_id' => $seq_id ) );
      $hijacked_segment = $self->_create_collection( \%args );
      $hijacked_segment->add_sequences( $sequence ) if defined( $sequence );
      $hijacked_segment->add_features(
        values %{ $self->{ '_seq_id_to_feature_table' }->
                         { $seq_id }->
                         { '_identifiable_features' }
                }
      );
      foreach my $start ( keys %{ $self->{ '_seq_id_to_feature_table' }->
                                         { $seq_id }->
                                         { '_anonymous_features' }
                                } ) {
        $hijacked_segment->add_features(
          @{ $self->{ '_seq_id_to_feature_table' }->
                    { $seq_id }->
                    { '_anonymous_features' }->
                    { $start }
           }
        );
      } # End foreach $start

      # Reuse %args
      %args = @_;

      # If the user named the sequence itself in their request, give it to 'em.
      # This entails removing from the request any args that name
      # the sequence.
      ## TODO: Dehackify
      my ( $unique_ids, $namespace, $names );
      if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
        ( $unique_ids, $namespace, $names ) =
          rearrange(
            [ [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
              [ qw( NAMESPACE NAME_SPACE CLASS ) ],
              [ qw( NAME NAMES DISPLAY_NAME DISPLAY_NAMES DISPLAYNAME DISPLAYNAMES ) ],
            ],
            @_
          );
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
      if( $unique_ids ) {
        my @to_be_removed;
        for( my $i = 0; $i < scalar( @$unique_ids ); $i++ ) {
          if( ( $unique_ids->[ $i ] eq $seq_id ) ||
              ( $namespace.':'.$unique_ids->[ $i ] eq $seq_id ) ||
              ( $unique_ids->[ $i ] eq $namespace.':'.$seq_id ) ) {
            ## TODO: REMOVE
            #warn "\$unique_ids->[ $i ] is $seq_id" if Bio::Graphics::Browser::DEBUG;
            push( @to_be_removed, $i );
          }
        }
        my $offset = 0;
        foreach my $i ( @to_be_removed ) {
          splice( @$unique_ids, $i + $offset++, 1 );
        }
        # If we removed any unique_ids, modify @_.
        if( $offset ) {
          if( scalar( @$unique_ids ) ) {
            $args{ '-unique_ids' } = $unique_ids;
          } else {
            delete $args{ '-unique_ids' };
          }
          # Oh yeah and make sure it's the only unique_ids argument..
          ## TODO: What if they're uppercase?
          # [ qw( UNIQUE_IDS UNIQUEIDS IDS UNIQUE_ID UNIQUEID ID ) ],
          delete @args{ qw( -unique_id -ids -id -uniqueids -uniqueid
                            uniqueids unique_id ids id uniqueids uniqueid ) };
          ## TODO: REMOVE
          #warn "\$unique_ids is now [ ".join( ', ', @$unique_ids )." ]" if Bio::Graphics::Browser::DEBUG;
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
      if( $names ) {
        my @to_be_removed;
        for( my $i = 0; $i < scalar( @$names ); $i++ ) {
          if( ( $names->[ $i ] eq $seq_id ) ||
              ( $namespace.':'.$names->[ $i ] eq $seq_id ) ||
              ( $names->[ $i ] eq $namespace.':'.$seq_id ) ) {
            ## TODO: REMOVE
            #warn "\$names->[ $i ] is $seq_id" if Bio::Graphics::Browser::DEBUG;
            push( @to_be_removed, $i );
          }
        }
        my $offset = 0;
        foreach my $i ( @to_be_removed ) {
          splice( @$names, $i + $offset++, 1 );
        }
        # If we removed any names, modify @_.
        if( $offset ) {
          if( scalar( @$names ) ) {
            $args{ '-names' } = $names;
          } else {
            delete $args{ '-names' };
          }
          # Oh yeah and make sure it's the only names argument..
          ## TODO: What if they're uppercase?
          delete @args{ qw( -name -display_name -display_names
                            -displaynames -displaynames
                            name names display_name display_names
                            displaynames displaynames ) };
          ## TODO: REMOVE
          #warn "\$names is now [ ".join( ', ', @$names )." ]" if Bio::Graphics::Browser::DEBUG;
        }
      }

      # TODO: Put this back, for efficiency.
      unless( %args ) {
        ## TODO: REMOVE
        #warn "hijacked segment is $hijacked_segment.  It has " . $hijacked_segment->feature_count() . " features.  We're returning it as-is." if Bio::Graphics::Browser::DEBUG;
        return $hijacked_segment;
      }

      ## TODO: REMOVE
      #warn "hijacked segment is $hijacked_segment.  It has " . $hijacked_segment->feature_count() . " features" if Bio::Graphics::Browser::DEBUG;
      #my $s = $hijacked_segment->get_collection( %args );
      #warn "Passing ( ". join( ', ', ( my @foo = %args ) ). " ) to its get_collection(..) method, we get " . ( ( defined $s ) ? $s : 'undef' ) . ", which has " . ( ( defined $s ) ? $s->feature_count() : 'N/A' ) . " features" if Bio::Graphics::Browser::DEBUG;
      #if( defined $s ) {
      #  my @sf = $s->features();
      #  warn "The first one of those is " . $sf[ 0 ] . ", which has " . $sf[ 0 ]->feature_count() . " subfeatures." if Bio::Graphics::Browser::DEBUG;
      #}
      #my @f = $hijacked_segment->features( %args );
      #warn "Passing ( ", join( ', ', ( my @foo = %args ) ), " ) to the hijacked segment's features(..) method, we get ( ", join( ', ', @f ), " )" if Bio::Graphics::Browser::DEBUG;
      #if( @f ) {
      #  warn "The first one of those has " . $f[ 0 ]->feature_count() . " features." if Bio::Graphics::Browser::DEBUG;
      #}
      #if( defined( $s ) ) {
      #  push( @segments, $s );
      #}

      # Now delegate to it & save the result, which will also be a Segment.
      push( @segments, $hijacked_segment->get_collection( %args ) );
    } # End foreach $seq_id

    if( wantarray ) {
      return @segments;
    } elsif( scalar( @segments ) > 1 ) {
      $self->throw( "Those features do not share a common sequence" );
    } else {
      return $segments[ 0 ];
    }
  }
} # get_collection(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 seq_ids

 Title   : seq_ids
 Usage   : my @seq_ids = $segmentprovider->seq_ids();
           OR
           my %seq_ids_and_counts =
               $segmentprovider->seq_ids( -count => 1 );
 Function: Enumerate all root seq_ids of features provided by this
           provider, and all seq_ids of sequences provided by this
           provider, and possibly count the features with each seq_id.
 Returns : a list of strings
           OR
           a hash mapping seq_id strings to integer counts
 Args    : see below

This routine returns a list of feature root seq_ids known to the
provider.  If the -count argument is given, it returns a hash of known
seq_ids mapped to their occurrence counts in this provider.  Note that
the returned list (or the keys of the returned hash) may include
seq_ids for which the count is 0, which indicates that the sequence is
provided but there are no features on it.

Arguments are -option=E<gt>value pairs as follows:

  -count aka -enumerate  if true, count the features

=cut

sub seq_ids {
  my $self = shift;
  my ( $count ) =
    rearrange( [ [ qw( COUNT ENUMERATE ) ] ], @_ );

  if( $count ) {
    my %feature_seq_ids =
      $self->Bio::DB::SimpleFeatureProvider::seq_ids( -count => 1 );
    # Add 0 entries for any unrepresented sequences.
    my @sequence_seq_ids = $self->Bio::DB::SimpleSequenceProvider::seq_ids();
    foreach my $seq_id ( @sequence_seq_ids ) {
      unless( defined( $feature_seq_ids{ $seq_id } ) ) {
        $feature_seq_ids{ $seq_id } = 0;
      }
    };
    return %feature_seq_ids;
  } else {
    my %intersection;
    my @feature_seq_ids = $self->Bio::DB::SimpleFeatureProvider::seq_ids();
    my @sequence_seq_ids = $self->Bio::DB::SimpleSequenceProvider::seq_ids();
    foreach my $seq_id ( @feature_seq_ids )  { $intersection{ $seq_id }++; };
    foreach my $seq_id ( @sequence_seq_ids ) { $intersection{ $seq_id }++; };
    return keys %intersection;
  }
} # seq_ids(..)

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
that is returned by this method creates collections of its own type
when its get_collection() method is called.  This is the likely
behavior anyway, but it should be noted, just in case.

=cut

# This implementation makes a new SimpleSegment with the given features.  It
# ignores the args except for '-seq_id' -- it will pass the seq_id to the new
# segment.
sub _create_collection {
  my $self = shift;
  my $args = shift;
  ## HACK because sometimes the args are passed in as a hash and
  ## sometimes as a list.
  if( $args && ( ref( $args ) eq 'ARRAY' ) ) {
    my %args_hash = @$args;
    $args = \%args_hash;
  }
  my @features = @_;
  if( @features && ( $features[ 0 ] eq 'lookup' ) ) {
    @features = $self->features( my @args = %$args );
  }
  ## TODO: REMOVE
  #warn "SimpleSegmentProvider::_create_collection( { ".join( ', ', ( my @foo = %$args ) )." }, ( ".join( ', ', @features )." )" if Bio::Graphics::Browser::DEBUG;

  ## TODO: REMOVE
  #if( @features ) {
  #  warn "creating collection with these features: ( ", join( ', ', @features ), " )" if Bio::Graphics::Browser::DEBUG;
  #}

  my $seq_id;
  ## TODO: REMOVE
  if( $args && ( ref( $args ) ne 'HASH' ) ) {
    use Data::Dumper;
    $self->throw( "Expected first arg to be a hash, but its a ".ref( $args ).": ".Dumper( $args )."." );
  }
  if( $args && %$args ) {
    $seq_id = $args->{ '-seq_id' };
  }
  return Bio::SeqFeature::SimpleSegment->new(
           '-features' => \@features,
           '-parent' => $self,
           ( defined $seq_id ? ( '-seq_id' => $seq_id ) : () )
         );
} # _create_collection(..)

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

  If a new feature was added as a result of this call, and that new
  feature's abs_range is a L<Bio::PrimarySeqI> object, then that
  sequence will be inserted into the sequence store if it is not
  already there, or updated if it is.

=cut

# Enveloped to add the new feature's sequence to the sequence store.
sub _insert_feature {
  my $self = shift;
  my $feature = shift;

  my $inserted = $self->SUPER::_insert_feature( $feature );
  if( defined $inserted ) {
    if( (
         my $abs_range = $inserted->abs_range()
        )->isa( 'Bio::PrimarySeqI' ) ) {
      ## TODO: REMOVE
      #warn "_insert_feature( $feature ): its ancestor is a real sequence.\n";
      $self->_insert_or_update_sequence( $abs_range );
    }
  }
  return $inserted;
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

  If this call is successful, and the returned feature's abs_range is
  a L<Bio::PrimarySeqI> object, then that sequence will be inserted
  into the sequence store if it is not already there, or updated if it
  is.

=cut

# Enveloped to add the new feature's sequence to the sequence store.
sub _insert_or_update_feature {
  my $self = shift;
  my $feature = shift;

  my $r_feature = $self->SUPER::_insert_or_update_feature( $feature );
  if( defined $r_feature ) {
    if( (
         my $abs_range = $r_feature->abs_range()
        )->isa( 'Bio::PrimarySeqI' ) ) {
      ## TODO: REMOVE`
      #warn "_insert_or_update_feature( $feature ): its ancestor is a real sequence.\n";
      $self->_insert_or_update_sequence( $abs_range );
    }
  }
  return $r_feature;
} # insert_or_update_feature(..)

1;

__END__
