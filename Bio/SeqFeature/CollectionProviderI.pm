# $Id$
#
# BioPerl module for Bio::SeqFeature::CollectionProviderI
#
# Cared for by Robert Hubley <rhubley@systemsbiology.org>
#
# Copyright Robert Hubley
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::CollectionProviderI - A provider of SeqFeatureI
objects as Bio::SeqFeature::CollectionIs that can update/add/remove
SeqFeatureIs passed in CollectionIs.

=head1 SYNOPSIS

    # get a CollectionProviderI somehow, perhaps a SimpleCollectionProvider
    my $provider = new Bio::SeqFeature::SimpleCollectionProvider();
    my $fg_color = $provider->get('fgcolor');

=head1 DESCRIPTION

A CollectionProviderI is an object that can return Bio::SeqFeature::CollectionIs of Bio::SeqFeatureIs, and can optionally update/add/remove SeqFeatures.

Features can be filtered by the following attributes:

  1) their range, perhaps relative to a "base" range, with a choice
     between features that are overlapping, contained within, or
     completely containing the given range

  2) their type

  3) other attributes using tag/value semantics

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Robert Hubley

Email rhubley@systemsbiology.org

=head1 CONTRIBUTORS

Paul Edlefsen, pedlefsen@systemsbiology.org
Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqFeature::CollectionProviderI;
use vars qw( @ISA );
use strict;
use Bio::Root::RootI;
use Carp;

@ISA = qw( Bio::Root::RootI );

#                   --Coders beware!--
# Changes to this get_collection() pod need to be copied to the following
# places, perhaps respecting existing modifications:
#   Bio/SeqFeature/SimpleCollectionProvider.pm [no  modifications]
#   Bio/DB/SegmentProviderI.pm                 [yes modifications]
# Also check up on
#   Bio/SeqFeature/CollectionI.pm
# , which has a similar interface for its features() method.

=head2 get_collection

 Title   : get_collection
 Usage   : my $collection = $collectionprovider->get_collection( %args );
           OR
           my $collection = $collectionprovider->get_collection( @types );
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
auto-sync (they might!).

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
                 default -baserange.  If this CollectionProviderI is also a
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
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 insert_or_update_collection

 Title   : insert_or_update_collection
 Usage   : $collectionprovider->insert_or_update($collection);
 Function: Attempts to update all the features of a collection.  If
           a feature doesn\'t exist it inserts it automatically.
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

sub insert_or_update_collection {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 insert_collection

 Title   : insert_collection
 Usage   : $collectionprovider->insert_collection($collection);
 Function: Insert all the features of a collection.  If any features
           already exist throw an exception. 
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

sub insert_collection {
   my ($self) = @_;
   $self->throw_not_implemented();
}

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
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 remove_collection

 Title   : remove_collection
 Usage   : $collectionprovider->remove_collection($collection);
 Function: Removes all the features in a collection.  If any features 
           do not exists throw an exception.
 Returns : None
 Args    : L<Bio::SeqFeature::CollectionI> object

=cut

sub remove_collection {
  my ($self) = @_;
  $self->throw_not_implemented();
}

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
  shift->throw_not_implemented;
}

=head2 parent_collection_provider

 Title   : parent_collection_provider
 Usage   : my $parent = $collectionprovider->parent_collection_provider();
 Function: Return the CollectionProviderI that is the parent of this provider.
 Returns : a L<Bio::SeqFeature::CollectionProviderI> or undef if there is none
 Args    : none

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
  shift->throw_not_implemented;
}

=head2 factory

 Title   : factory
 Usage   : my $parent = $collectionprovider->factory();
 Function: Return the CollectionProviderI that is the parent of this provider.
 Returns : a L<Bio::SeqFeature::CollectionProviderI> or undef if there is none
 Args    : none

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

  This method is implemented in the interface as a (glob ref) alias to
  parent_collection_provider().

=cut

  *factory = \&parent_collection_provider;

1;

__END__
