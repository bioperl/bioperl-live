package Bio::DB::SegmentProviderI;

# $Id $
# An interface for objects that are both Bio::DB::SequenceProviderIs and
# Bio::DB::FeatureProviderIs.

=head1 NAME

Bio::DB::SegmentProviderI -- A provider of collections of sequence
features from a database or other non-trivial backing store that also
provides the sequences that the features reside on.

=head1 SYNOPSIS

## TODO: Update this for Segments (it's copied from FeatureProvider).

 use Bio::DB::SimpleSegmentProvider;

 use Bio::SeqFeature::Generic;
 use Bio::SeqFeature::SimpleCollection;

 my $data_provider =
   Bio::SeqFeature::SimpleSegmentProvider->new();

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
Bio::SeqFeature::SegmentIs stored in a database or other (generally
external) backing store.  It is a L<Bio::DB::FeatureProviderI> and a
L<Bio::DB::SequenceProviderI>.

The only new method is segment(), which is an alias for
get_collection().  The get_collection() method of a SegmentProviderI
returns a L<Bio::SeqFeature::SegmentI> object (which is an additional
constraint over the requirement from the L<Bio::DB::FeatureProviderI>
interface, which requires that the get_collection() method return any
L<Bio::SeqFeature::CollectionI> object).

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

use Bio::DB::FeatureProviderI;
use Bio::DB::SequenceProviderI;
@ISA = qw( Bio::DB::FeatureProviderI Bio::DB::SequenceProviderI );

use vars '$VERSION';
$VERSION = '1.00';

#                   --Coders beware!--
# This pod is a modification of the pod for get_collection() in
#   Bio/SeqFeature/CollectionProviderI.pm
# , so changes must be kept in sync.
# Also note that
#   Bio/SeqFeature/SegmentI.pm
# copies and modifies this pod, changes must be kept in sync.

=head2 get_collection

 Title    : get_collection
 Usage    : my $segment = $segmentprovider->get_collection( %args );
            OR
            my $segment = $segmentprovider->get_collection( @types );
 Returns  : A L<Bio::SeqFeature::SegmentI> object
 Args     : see below
 Status   : Public
 Exception: "Those features do not share a common sequence" if the
            features that would otherwise be included in the resulting
            segment do not all fall on the same sequence.

NOTE: This method is (almost) identical to the get_collection() method
from L<Bio::DB::FeatureProviderI> that it overrides.  The entire
documentation follows, but first a brief summary of the changes:
  * This method returns L<Bio::SeqFeature::SegmentI> objects instead
    of mere CollectionI objects.  SegmentI objects are CollectionI
    objects, so this is an additional constraint on the interface.
    The returned SegmentI objects will have as their range the range
    searched, if any, or the smallest range that encloses the returned
    features.
  * This method will throw an exception if the features that would be
    included in the resulting SegmentI do not all share a common
    sequence. The returned SegmentI will have as its seq_id the common
    sequences' unique_id() or primary_id().

This routine will retrieve a L<Bio::SeqFeature::SegmentI> object based
on feature type, location or attributes.  The SeqFeatureI objects in
the returned SegmentI may or may not be newly instantiated by this
request.  They will have as their range the range searched, if any, or
the smallest range that encloses the returned features.  They will
have as their seq_id() the unique_id() or primary_id() of the returned
features' common sequence.  If the returned features do not share a
common sequence then an exception will be thrown.

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

sub get_collection {
  shift->throw_not_implemented();
}

=head2 segment

 Title   : segment
 Usage   : my $segment = $segmentprovider->segment( %args );
           OR
           my $segment = $segmentprovider->segment( @types );
 Returns : A L<Bio::SeqFeature::SegmentI> object
 Args    : see below
 Status  : Public

  This method is a (glob ref) alias for get_collection().

=cut

  *segment = \&get_collection;

1;

__END__
