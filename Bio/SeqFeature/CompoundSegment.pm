package Bio::SeqFeature::CompoundSegment;

# $Id$
# An implementation of Bio::SeqFeature::SegmentI that provides the
# data (the features and sequences) from another implementation in
# addition to its own.

=head1 NAME

Bio::SeqFeature::CompoundSegment -- An implementation of
Bio::SeqFeature::SegmentI that provides the data (the features and
sequences) from another implementation in addition to its own.

=head1 SYNOPSIS

TODO

=head1 DESCRIPTION

TODO

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

use Bio::SeqFeature::SimpleSegment;
use Bio::DB::CompoundSegmentProvider;
@ISA = qw( Bio::SeqFeature::SimpleSegment
           Bio::DB::CompoundSegmentProvider );

use vars '$VERSION';
$VERSION = '1.00';

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'
use Bio::RelRange qw( &absSeqId );
use Bio::SeqFeature::CompoundIterator;

sub new {
  my( $class, @args ) = @_;

  my $self = $class->SUPER::new( @args );
  $self->_initialize_compound_segment( @args );
  return $self;
} # new(..)

sub _initialize_compound_segment {
  my $self = shift;
  my @args = @_;

  return if( $self->{ '_compound_segment_initialized' } );

  $self->_initialize_compound_segment_provider( @args );
  $self->_initialize_simple_segment( @args );

  $self->{ '_compound_segment_initialized' }++;
  return $self;
} # _initialize_compound_segment(..)

## TODO: Document this.
sub add_next_segment {
  my $self = shift;

  foreach my $next_segment ( @_ ) {
    ## TODO: REMOVE
    #warn "Adding next segment: $next_segment";
    push( @{ $self->{ '_next_providers' } }, $next_segment );
  }
  return;
} # add_next_segment(..)

## TODO: Document this.
sub get_next_segments {
  my $self = shift;

  return grep { $_->isa( 'Bio::SeqFeature::SegmentI' ) }
              @{ $self->{ '_next_providers' } };
} # get_next_segments()

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

# Delegates to Bio::DB::CompoundSegmentProvider.
sub get_collection {
  my $self = shift;

  return $self->Bio::DB::CompoundSegmentProvider::get_collection( @_ );
} # get_collection(..)

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
  my ( $iterator, $callback );
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $iterator, $callback ) =
      rearrange(
        [ [ qw( ITERATOR STREAM ) ],
          [ qw( CALLBACK CALL_BACK ) ]
        ],
        @_
      );
  }
  # There's 3 options for this: by default features() returns a
  # list. if -iterator is given then it returns an iterator.  if -callback
  # is given then some method is called.  We need to deal with each option
  # separately.
  if( $iterator ) {
    my $compound_iterator = Bio::SeqFeature::CompoundIterator->new();
    $compound_iterator->add_next_iterator( $self->SUPER::features( @_ ) );
    foreach my $next_segment ( $self->get_next_segments() ) {
      $compound_iterator->add_next_iterator( $next_segment->features( @_ ) );
    }
    return $compound_iterator;
  } elsif( $callback ) {
    # Callbacks return false when they fail, and once they fail the
    # callbacks should cease.
    if( $self->SUPER::features( @_ ) ) {
      foreach my $next_segment ( $self->get_next_segments() ) {
        last unless $next_segment->features( @_ );
      }
    }
    return 1;
  } else {
    ## TODO: Put back.  Testing.
    my @features;# = $self->SUPER::features( @_ );
    foreach my $next_segment ( $self->get_next_segments() ) {
      ## TODO: REMOVE
      warn "CompoundSegment: next_segment is $next_segment, a ".ref( $next_segment )."\n";
      push( @features, $next_segment->features( @_ ) );
    }
    return @features;
  } # End if $iterator .. or $callback .. or else.
} # features(..)

sub _non_compound_features {
  return shift->SUPER::features( @_ );
} # _non_compound_features(..)

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
  my $count = $self->SUPER::feature_count( @_ );
  foreach my $next_segment ( $self->get_next_segments() ) {
    $count += $next_segment->feature_count( @_ );
  }
  return $count;
} # feature_count()

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 insert_or_update_collection

 Title   : insert_or_update_collection
 Usage   : $collectionprovider->insert_or_update($collection);
 Function: Attempts to update all the features of a collection.  If
           a feature doesn\'t exist it inserts it automatically.
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

sub insert_or_update_collection {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::insert_or_update_collection( @_ );
} # insert_or_update_collection(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

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

  $self->Bio::DB::CompoundSegmentProvider::insert_collection( @_ );
} # insert_collection(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

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

  $self->Bio::DB::CompoundSegmentProvider::update_collection( @_ );
} # update_collection(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 remove_collection

 Title   : remove_collection
 Usage   : $collectionprovider->remove_collection($collection);
 Function: Removes all the features in a collection.  If any features 
           do not exists throw an exception.
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

sub remove_collection {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::remove_collection( @_ );
} # remove_collection(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

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

=cut

sub types {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::types( @_ );
} # types(..)

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

  $self->Bio::DB::CompoundSegmentProvider::seq_ids( @_ );
} # seq_ids(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 get_all_primary_ids

 Title   : get_all_primary_ids
 Usage   : my @primary_ids = $provider->get_all_primary_ids()
 Function: Returns an array of all the primary_ids of the sequence
           objects in this data store. These maybe ids (display style)
           or accession numbers or something else completely different
           - they may be anything so long as each sequence has a
           different one.  Note that although some sequences may have
           undefined primary_ids (bad!), the returned list will not
           include undef.
 Returns : an array of strings
 Args    : none
 Status  : Public

=cut

sub get_all_primary_ids {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::get_all_primary_ids( @_ );
} # get_all_primary_ids(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 unique_ids

 Title     : unique_ids
 Usage     : my @unique_ids = $sequenceprovider->unique_ids();
 Function  : Return a list of the unique_ids of all sequences provided
             by this SequenceProvider.  Note that although some
             sequences may have undefined unique_ids, the returned
             list will not include undef.
 Returns   : an array of strings.
 Args      : none
 Status    : Public

=cut

sub unique_ids {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::unique_ids( @_ );
} # unique_ids(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 sequence_count

 Title   : sequence_count
 Usage   : $collection->sequence_count()
 Function: Return the number of L<Bio::PrimarySeqI> objects that would
           be returned by a call to sequences() with no arguments.
 Returns : integer representing the number of sequence objects
 Args    : None

  This method is implemented in the interface to return
    scalar( $self->sequences() )
  Because this is not particularly efficient, implementers are
  encouraged to override it, but the result should of course be the
  same.

=cut

sub sequence_count {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::sequence_count( @_ );
} # sequence_count()

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 add_sequences

 Title   : add_sequences
 Usage   : my @added = $collection->add_sequences( @sequence_list );
 Function: Adds the given sequences to this provider.
 Returns : The sequences added (or their count, in scalar context).
 Args    : An array of L<Bio::PrimarySeqI>s
 Status  : Public

=cut

sub add_sequences {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::add_sequences( @_ );
} # add_sequences(..)

#                   --Coders beware!--
# Changes to the interface pod need to be copied to here.

=head2 remove_sequences

 Title   : remove_sequences
 Usage   : my @removed = $collection->remove_sequences( @sequence_list )
 Function: Removes the requested sequences from this provider.
 Returns : The removed sequences (or their count, in scalar context)
 Args    : An array of L<Bio::PrimarySeqI>s or their ids (see below)
 Status  : Public

  If any argument is a string, it will be taken as either the
  unique_id or the primary_id or the accession of a sequence to be
  removed.  The return list will contain the actual removed sequence
  if there was a match.

=cut

sub remove_sequences {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::remove_sequences( @_ );
} # remove_sequences()

=head2 sequences

 Title   : sequences
 Usage   : my @seqs = $provider->sequences( @names );
           OR
           my @seqs = $provider->sequences( %args );
 Function: Retrieves a list of L<Bio::PrimarySeqI> objects.
 Returns : a list of L<Bio::PrimarySeqI> objects
           OR
           (when the -iterator option is true)
             a L<Bio::Seq::IteratorI> object
           OR
           (when the -callback option is true) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public
 Throws  : "$some_text does not exist" ($some_text might be anything)
            if a version is given and a matching sequence exists,
            but not of that version.

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of sequence names (as if they
were given as -names => \@_).  In the named parameter form, the
arguments are a series of -name=E<gt>value pairs.  Note that the table
below is not exhaustive; implementations must support these but may
support next arguments as well (and are responsible for documenting the
difference).

  Argument       Description
  --------       ------------

  -unique_id     A (string) unique_id.  See also -namespace and -id.
  -unique_ids    An array reference to multiple unique_id values.

  -primary_id    A (string) primary_id.  See also -namespace and -id.
  -primary_ids   An array reference to multiple primary_id values.

  -display_id    A (string) display_id.  See also -namespace and -id.
  -display_ids   An array reference to multiple display_id values.

  -id            A (string) unique_id or primary_id or display_id.
                 See also -namespace.
  -ids           An array reference to multiple id values.

  -accession     A (string) accession.  See also -namespace and -name.
  -accessions    An array reference to multiple accession values.

  -name          A (string) accession or id.  See also -namespace and -id.
  -names         An array reference to multiple accession/id values.

  -namespace     A (string) namespace qualifier to help resolve the names.
  -class         same as -namespace

  -iterator      Return a L<Bio::Seq::IteratorI>
  -stream        same as -iterator

  -callback      A callback to invoke on each sequence

All plural arguments are interchangeable with their singular counterparts.

The -unique_ids argument is a reference to a list of strings.  Every
returned sequence must have its unique_id value in this list if the
list is provided*.  Iff a sequence has no unique_id value then its
primary_id will be used instead.  A -unique_id argument is treated as
a single-element list of unique_ids.

The -ids argument is a reference to a list of strings.  Every returned
sequence must have its unique_id, primary_id, or display_id value in
this list if the list is provided*.  An -id argument is treated as a
single-element list of ids.

The -accessions argument is a reference to a list of strings.  Every
returned sequence must have its accession value in this list if the
list is provided*.  An -accession argument is treated as a
single-element list of names.  If the accession value contains a dot
('.') then the suffix after the final dot will be interpreted as a
version number.  If the accession is available but not in the given
version then an exception (ending in ' does not exist') will be
thrown.  An empty version value is acceptable and means 'latest version',
which is also the default, so -accession => 'foo' and -accession =>
'foo.' mean the same thing.

NOTE: If your accession value contains a dot (unrelated to version
number) then you should postpend the value with a dot to disambiguate:
'accession.that.contains.dots' will be interpreted as the accession
'accession.that.contains' with version 'dots' unless you postpend it
with a dot, like so: -accession => 'accession.that.contains.dots.'

The -names argument is a reference to a list of strings.  Every
returned sequence must have its accession, unique_id, primary_id,
or display_id value in this list if the list is provided*.  A -name
argument is treated as a single-element list of names.

If a -namespace is provided then names and ids (both queries and
targets) will be prepended with "$namespace:" as a bonus.  So if you
do sequences( -names => [ 'foo', 'bar' ], -namespace => 'ns' ) then any
sequence with the accession, unique_id, primary_id, or display_id
'foo', 'ns:foo', 'bar', or 'ns:bar' will be returned.

If -iterator is true, then the method returns an object of type
L<Bio::Seq::IteratorI>.  Each call to next_seq() on this
object returns a L<Bio::PrimarySeqI> object from this provider.

If -callback is passed a code reference, the code reference will be
invoked on each sequence returned.  The code will be passed two
arguments consisting of the current sequence and this SequenceProviderI
object, and must return a true value.  If the code returns a false
value, sequence retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

Footnote (*): -unique_ids, -primary_ids, -display_ids, -ids,
-accessions, and -names will be ORed together if they appear together,
so if you do $provider->sequences( -unique_id => 'foo', -accession =>
'bar' ) then you will get all sequences with unique_id 'foo' and
(also) all sequences with accession 'bar'.

NOTE: the following methods all build on top of sequences(), and do not
need to be explicitly implemented.

    get_Seq_by_id()
    get_Seq_by_primary_id()
    get_Seq_by_accession()
    get_Seq_by_version()
    get_Seq_stream()
    get_PrimarySeq_stream()

=cut

sub sequences {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::sequences( @_ );
} # sequences(..)

=head2 sequence_count

 Title   : sequence_count
 Usage   : $collection->sequence_count()
 Function: Return the number of L<Bio::PrimarySeqI> objects that would
           be returned by a call to sequences() with no arguments.
 Returns : integer representing the number of sequence objects
 Args    : None

=cut

sub sequence_count {
  my $self = shift;

  $self->Bio::DB::CompoundSegmentProvider::sequence_count( @_ );
} # sequence_count()

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
  shift->_create_segment( @_ );
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
  my ( $args, @features ) = @_;
  ## HACK because sometimes the args are passed in as a hash and
  ## sometimes as a list.
  if( $args && ( ref( $args ) eq 'HASH' ) ) {
    my @args_list = %$args;
    $args = \@args_list;
  }
  if( @features && ( $features[ 0 ] eq 'lookup' ) ) {
    @features = $self->_non_compound_features( @$args );
  }

  ## TODO: REMOVE
  #warn "CompoundSegment::_create_segment( { ".join( ', ', @$args )." }, ( ".join( ', ', @features )." )" if Bio::Graphics::Browser::DEBUG;

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
    #warn "The baserange that was given, $baserange, isa ".ref( $baserange ).", but we need a RangeI.  Using \$self->abs_range().";
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
  #warn "CompoundSegment::_create_segment(..): union_range is $union_range, baserange is " . $baserange->Bio::RelRangeI::toString() . '.' if Bio::Graphics::Browser::DEBUG;
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
    my $abs_baserange = $baserange->abs_range();
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
  #warn "CompoundSegment::_create_segment(..): \%args_to_new are { ".join( ', ', ( my @foo = %args_to_new ) )." }" if Bio::Graphics::Browser::DEBUG;
  my $new_segment = $self->new( %args_to_new );

  ## TODO: Add _create_collection to the interface def for CollectionProviderI.
  foreach my $next_provider ( $self->get_next_providers() ) {
    $new_segment->add_next_provider(
      $next_provider->_create_collection( @_ )
    );
  }
  return $new_segment;
} # _create_segment(..)

## TODO: REMOVE.  Testing.
sub toString {
  return shift->Bio::SeqFeature::SimpleSegment::toString( @_ );
}

1;

__END__
