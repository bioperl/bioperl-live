package Bio::SeqFeature::SimpleCollection;

# $Id $

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

=head2 add_features

 Title   : add_features
 Usage   : $collection->add_features( @feature_list );
 Function: Adds the given features to this Collection.
 Returns : The features added (or their count, in scalar context).
 Args    : An array of SeqFeatures or their ids
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
 Args    : An array of SeqFeatures or their ids
 Status  : Public

=cut

sub remove_features {
    shift->throw_not_implemented();
}

=head2 features

 Title   : features
 Usage   : @features = $collection->features( @args );
 Returns : a list of Bio::SeqFeatureI objects, or an iterator
 Args    : see below
 Status  : Public

This routine will retrieve features associated with this collection
object.  It can be used to return all features, or a subset based on
their type, location, or attributes.

If the implementing class also implements Bio::RangeI, then the feature
locations will use coordinates relative to the reference sequence in
effect at the time that features() was called.

If a range is specified using either the -range or the -start, -end,
and -strand arguments, then this range will be used to narrow the
results, according to the specified -rangetype and -strandtype
arguments.

-rangetype is one of:
   "overlaps"      return all features that overlap the range (default)
   "contains"      return features completely contained within the range
   "contained_in"  return features that completely contain the range

Note that if the implementing class implements RangeI then the range
will default to that range.  If it implements RangeI and a range is
specified as an argument, then the coordinates of the given range will
be interpreted relative to the implementing class\'s range.  If the
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

  -iterator      Return a Bio::SeqFeature::IteratorI

  -callback      A callback to invoke on each feature

The -attributes argument is a hashref containing one or more
attributes to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple string matching, and multiple attributes
are ANDed together.  More complex filtering can be performed using the
-callback option (see below).

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

NOTE: the following methods all build on top of features(), and do not
need to be explicitly implemented.

    features_in_range()
    overlapping_features()
    contained_features()
    contained_in()
    get_feature_stream()

=cut

sub features {
  my $self = shift;
  my @features_to_return;
  ## TODO: Something more sophisticated.  This always returns all features.
  push( @features_to_return, values $self->{ '_identifiable_features' } );
  foreach my $start ( keys $self->{ '_anonymous_features' } ) {
    push( @features_to_return,
          @{ $self->{ '_anonymous_features' }{ $start } }
        );
  }
  return @features_to_return;
} # features()

1;

__END__
