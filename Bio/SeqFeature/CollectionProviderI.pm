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

=head2 get_collection

 Title   : get_collection
 Usage   : my $collection = $collectionprovider->get_collection(@args);
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

If a range is specified using the -range argument then this range will
 be used to narrow the results, according to the specified -rangetype
 and -strandtype arguments.

-rangetype is one of:
   "overlaps"      return all features that overlap the range (default)
   "contains"      return features completely contained within the range
   "contained_in"  return features that completely contain the range

Note that if the implementing class implements RangeI then the baselocation
will default to that range.  If a baselocation is given or defaulted and a range is specified as an argument, then the coordinates of the given range will
be interpreted relative to the implementing class's range.  If the '
implementing class does not implement RangeI and no range is given,
then -rangetype may be ignored.

-strandmatch is one of:
   "strong"        ranges must have the same strand
                   (default ONLY when -strand is specified and non-zero)
   "weak"          ranges must have the same strand or no strand
   "ignore"        ignore strand information
                   (default unless -strand is specified and non-zero)

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of feature types.  In the
named parameter form, the arguments are a series of -name=E<gt>value
pairs.

  Argument       Description
  --------       ------------

  -type          A type name or an object of type Bio::SeqFeature::TypeI
  -types         An array reference to multiple type names or TypeI objects

  -attributes    A hashref containing a set of attributes to match.  See
                 below.

  -location      A Bio::LocationI object defining the range to search and
                 the rangetype.  Use -range (and -baselocation,
                 perhaps; see below) as an alternative to -location.
                 See also -strandmatch.  There may be a default value
                 for -location.

  -baselocation  A Bio::LocationI object defining the location to which
                 the -range argument is relative.  There may be a
                 default -baselocation.  If this CollectionI is also a
                 Bio::RangeI, then the default -baselocation should be
                 its range.

  -range         A Bio::RangeI object defining the range to search.  See also
                 -strandmatch and -rangetype.  Use instead of
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

=cut

sub get_collection {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 insert_or_update_collection

 Title   : insert_or_update_collection
 Usage   : $collectionprovider->insert_or_update($collection);
 Function: Attempts to update all the features of a collection.  If
           a feature doesn't exist it inserts it automatically. '
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

1;
