# $Id$
#
# BioPerl module for Bio::SeqFeature::CollectionI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::CollectionI - An interface for a collection of
SeqFeatureI objects.

=head1 SYNOPSIS


# get a Bio::SeqFeature::CollectionI somehow
# perhaps a Bio::SeqFeature::SimpleCollection

 use Bio::SeqFeature::SimpleCollection;
 my $collection = new Bio::SeqFeature::SimpleCollection;
 $collection->add_features( @featurelist );
 
 $collection->features(-range => new Bio::Location::Simple
                                        (-start=> 1, -end => 300),
                       -rangetype=>'overlaps');

=head1 DESCRIPTION

This interface describes the basic methods needed for a collection of Sequence Features.  

Features can be filtered by the following attributes:

  1) their range, perhaps relative to a "base" range, with a choice
     between features that are overlapping, contained within, or
     completely containing the given range

  2) their type

  3) other attributes using tag/value semantics

Access to the contained features can be achieved using the
CollectionProviderI interface to produce a new semantic view on this
collection, or via accessing the feature list directly.  Access to the
feature list can be achieved using multiple techniques:

  1) as another CollectionI

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::CollectionI;
use vars qw( @ISA );
use strict;
use Carp;
use Bio::Root::RootI;
use Bio::FeatureHolderI;;
use Bio::SeqFeature::CollectionProviderI;

@ISA = qw( Bio::Root::RootI
           Bio::FeatureHolderI
           Bio::SeqFeature::CollectionProviderI
         );


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
what it was beforehand.

=cut

sub sorted {
  shift->throw_not_implemented();
}

=head2 add_features

 Title   : add_features
 Usage   : $collection->add_features( @feature_list );
 Function: Adds the given features to this Collection.
 Returns : The features added (or their count, in scalar context).
 Args    : An array of L<Bio::SeqFeatureI>s
 Status  : Public

=cut

sub add_features {
    shift->throw_not_implemented();
}

=head2 remove_features

 Title   : remove_features
 Usage   : $collection->remove_features( @feature_list )
 Function: Removes the requested sequence features
 Returns : The removed features (or their count, in scalar context)
 Args    : An array of L<Bio::SeqFeatureI>s or their unique_ids
 Status  : Public

=cut

sub remove_features {
    shift->throw_not_implemented();
}

#                   --Coders beware!--
# Changes to this features() pod need to be copied to the following
# places, perhaps respecting existing modifications:
#   Bio/SeqFeature/SimpleCollection.pm [no  modifications]
#   Bio/SeqFeatureI.pm                 [no  modifications]
#   Bio/SeqFeature/SegmentI.pm         [yes modifications]
# Also check up on
#   Bio/SeqFeature/CollectionProviderI.pm
# , which has a similar interface for its get_collection() method.

=head2 features

 Title   : features
 Usage   : @features = $collection->features( %args );
           OR
           @features = $collection->features( @types );
 Returns : a list of L<Bio::SeqFeatureI> objects,
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
object returns a Bio::SeqFeatureI object from this collection.

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

sub features { shift->throw_not_implemented }

=head2 overlapping_features

 Title   : overlapping_features
 Usage   : @features = $collection->overlapping_features( %args )
 Function: get features that overlap the range of this collection
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : same as features()
 Status  : Public

This method is identical to features() except that it defaults to
finding overlapping features.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub overlapping_features {
  my $self = shift;
  if( $_[ 0 ] =~ /^-/ ) {
    return $self->features( @_, -rangetype => 'overlaps' );
  } else {
    return $self->features( -types => \@_, -rangetype => 'overlaps' );
  }
} # overlapping_features()

=head2 contained_features

 Title   : contained_features
 Usage   : @features = $collection->contained_features( %args )
 Function: get features that are contained in the range of this collection
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : same as features()
 Status  : Public

This method is identical to features() except that it defaults to
a range type of 'contains'.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub contained_features {
  my $self = shift;
  if( $_[ 0 ] =~ /^-/ ) {
    return $self->features( @_, -rangetype => 'contains' );
  } else {
    return $self->features( -types => \@_, -rangetype => 'contains' );
  }
} # contained_features()

=head2 contained_in

 Title   : contained_in
 Usage   : @features = $collection->contained_in( %args )
 Function: get features that contain the range of this collection
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : same as features()
 Status  : Public

This method is identical to features() except that it defaults to
a range type of 'contained_in'.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub contained_in {
  my $self = shift;
  if( $_[ 0 ] =~ /^-/ ) {
    return $self->features( @_, -rangetype => 'contained_in' );
  } else {
    return $self->features( -types => \@_, -rangetype => 'contained_in' );
  }
} # contained_in()

=head2 get_feature_stream

 Title   : get_feature_stream
 Usage   : $iterator = $collection->get_feature_stream( %args )
 Function: get an iterator over the features in this collection
 Returns : a Bio::SeqFeature::IteratorI
 Args    : same as features()
 Status  : Public

This method is identical to features() except that it always generates
an iterator.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub get_feature_stream {
  my $self = shift;
  if( $_[ 0 ] =~ /^-/ ) {
    return $self->features( @_, -iterator => 1 );
  } else {
    return $self->features( -types => \@_, -iterator => 1 );
  }
} # get_feature_stream()

=head2 features_in_range

 Title   : features_in_range
 Usage   : @features = $collection->features_in_range( $range );
             OR
           @features = $collection->features_in_range( %args );
 Function: Retrieves a list of features which were contained or overlap the
           the requested range
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : same as features(), or a single L<Bio::RangeI> argument
 Status  : Public

This method is identical to features() except that its first argument, if it is a RangeI object, will be used as the -range argument.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub features_in_range {
  my $self = shift;
  my @args = @_;
  ## eegads.  TODO: Test this!
  if( ( $_[ 0 ] !~ /^-/ ) && ( UNIVERSAL::isa( $_[ 0 ], 'Bio::RangeI' ) ) ) {
    my $range = shift;
    if( $_[ 0 ] =~ /^-/ ) {
      @args = ( @_, -range => $range );
    } else {
      @args = ( -types => \@_, -range => $range );
    }
  }
  return $self->features( @args );
} # features_in_range()

=head2 get_feature_by_name

 Title   : get_feature_by_name
 Usage   : my @features = $collection->get_feature_by_name( $name )
           OR
           my @features = $collection->get_feature_by_name( $namespace, $name )
           OR
           my @features = $collection->get_feature_by_name( %args )
 Function: fetch features by their name
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : the string name or array ref of names
           OR
           the string namespace and the string name or array ref of names
           OR
           a hash, same as features()
 Status  : Public

This method is identical to features() except that it can take as
unnamed arguments the -name/-names value (one argument) OR the
-namespace value and the -name/-names value (two arguments).

Again, here's the deal:
  1) one argument: the argument is treated as the -name/-names value
  2) two arguments: the arguments are treated as the -namespace value and
       the -name/-names value, in that order.
     (note: this uses _rearrange() so the first argument must not
     begin with a hyphen or it will be interpreted as a named
     argument).
  3) an args hash as with features()

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub get_feature_by_name {
  my $self = shift;
  my @args = @_;
  if( $_[ 0 ] !~ /^-/ ) {
    if( scalar( @args ) == 2 ) {
      my ( $namespace, $name ) = @_;
      @args = ( -namespace => $namespace, -name => $name );
    } elsif( scalar( @args ) == 1 ) {
      my $name = shift;
      @args = ( -name => $name );
    } else {
      $self->throw( "Illegal number of arguments to get_feature_by_name(..): "
                    . scalar( @args ) );
    }
  }
  return $self->features( @args );
} # get_feature_by_name(..)

=head2 get_feature_by_id

 Title   : get_feature_by_id
 Usage   : my @features = $collection->get_feature_by_id( $unique_id )
           OR
           my @features = $collection->get_feature_by_id( %args )
 Function: fetch features by their unique_ids
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : the string unique_id or array ref of unique_ids
           OR
           a list of string unique_ids
           OR
           a hash, same as features()
 Status  : Public

This method is identical to features() except that it can take as
an unnamed argument(s) the -unique_id/-unique_ids value(s).

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub get_feature_by_id {
  my $self = shift;
  my @args = @_;
  if( $_[ 0 ] !~ /^-/ ) {
    if( scalar( @args ) == 1 ) {
      my $id = shift;
      @args = ( -unique_ids => $id );
    } else {
      @args = ( -unique_ids => \@_ );
    }
  }
  return $self->features( @args );
} # get_feature_by_id(..)

=head2 get_feature_by_attribute

 Title   : get_feature_by_attribute
 Usage   : my @features = $collection->get_feature_by_attribute( %attrs )
           OR
           my @features = $collection->get_feature_by_attribute( $attrs_ref, %args )
 Function: fetch features by their attributes
 Returns : a list of L<Bio::SeqFeatureI> objects,
           OR
           an iterator (when the -iterator option is true)
           OR
           (when the -callback option is given) true iff the callbacks
             completed.
 Args    : a hash, as would be passed to features() as -attributes => \%attrs
           OR
           a hash ref, as would be passed to features() as -attributes => $ref,
             then some other args to the features() method.
 Status  : Public

This method is identical to features() except that it assumes that the
given hash is the value meant for the -attributes argument, or (if the
first argument does not begin with '-') that the first argument is
meant as the -attributes value and the rest of them are usual
features() arguments.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub get_feature_by_attribute {
  my $self = shift;
  my @args = @_;
  if( $_[ 0 ] !~ /^-/ ) {
    my $attributes = shift;
    @args = ( -attributes => $attributes, @_ );
  } else {
    @args = ( -attributes => \@_ );
  }
  return $self->features( @args );
} # get_feature_by_attribute(..)

1;

__END__
