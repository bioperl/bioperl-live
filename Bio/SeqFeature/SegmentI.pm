package Bio::SeqFeature::SegmentI;

# $Id $
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
@ISA = qw( Bio::RelRangeI
           Bio::SeqFeature::CollectionI
           Bio::DB::SegmentProviderI );

use vars '$VERSION';
$VERSION = '1.00';

#                   --Coders beware!--
# This pod is a modification of the pod for features() in
#   Bio/SeqFeature/CollectionI.pm
# , so changes must be kept in sync.

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
    coordinates that are relative.  Features that are returned in
    absolute mode will be returned with absolute coordinates.  The
    mode is determined by the -baserange and -absolute arguments and
    by the absolute() flag, in that precedence order.

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

If a -baserange is specified then unqualified ranges given with the
-ranges argument will be interpreted as relative to that baserange.
If no -baserange is given then a default will be provided that will
depend on the value of the absolute() flag.  If absolute() is true
then the default is the value returned by the abs_seq_id() method; if
absolute() is false then the default is this SegmentI object ($self).
Note that the baserange only applies to unqualified ranges, ie. ranges
that have no defined seq_id.  You may force absolute range
interpretations by giving a -baserange that is not a L<Bio::RangeI>
(such as the string 'absolute', though any string will do the trick),
by providing a true value to the -absolute argument, by setting the
absolute() flag to true, or by qualifying all given ranges.

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
                 baserange affects the sort order.

  -absolute      If -absolute is given and true then all behavior will be as
                 if this SegmentI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SegmentI's
                 absolute() flag was set to false, even if it isn't.
                 Note that -baserange can still be given and can force
                 absoluteness or relativeness, and that takes
                 precedence over -absolute.

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
  shift->throw_not_implemented();
}

#                   --Coders beware!--
# This pod is a modification of the pod for get_collection() in
#   Bio/DB/SegmentProviderI.pm
# , so changes must be kept in sync.
# Changes to this features() pod need to be copied to the following
# places, perhaps respecting existing modifications:
#   Bio/SeqFeature/SimpleSegment.pm [no  modifications]

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
                 baserange affects the sort order.

  -absolute      If -absolute is given and true then all behavior will be as
                 if this SegmentI's absolute() flag was set to true,
                 even if it isn't.  If -absolute is given and false
                 then all behavior will be as if this SegmentI's
                 absolute() flag was set to false, even if it isn't.
                 Note that -baserange can still be given and can force
                 absoluteness or relativeness, and that takes
                 precedence over -absolute.

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
  shift->throw_not_implemented();
}

1;

__END__
