package Bio::SeqFeature::AggregatorI.pm

# $Id $
# Factory for grouping features into components of a containing feature.

=head1 NAME

Bio::SeqFeature::AggregatorI -- Factory for grouping features into
components of a containing (parent) feature.

=head1 SYNOPSIS

  ## Get an AggregatorI somehow.  Perhaps Bio::DB::GFF::Aggregator::orf.
  use Bio::DB::GFF::Aggregator::orf;
  my $aggregator = new Bio::DB::GFF::Aggregator::orf();

  # If you have a list of features you may aggregate them like so:
  my $did_aggregate = $aggregator->aggregate( $features_array_ref );
  # The return value is true iff the aggregator modified the given list.
  # Note that if aggregation occurred, the list will have additional
  # features but the original features will not be removed.

  # Sometimes aggregators need to know the CollectionProviderI that is
  # responsible for providing the features.  You may specify it.
  $aggregator->aggregate( $features_array_ref, $provider );
  # Or you may give it a CollectionI object:
  my $did_aggregate_collection = $aggregator->aggregate( $collection );
  # The return value in this case will be true iff the given
  # collection was modified.

  # If you would like to know what types of features an aggregator
  # needs in order to do its aggregation, you may use
  # the disaggregate_types() method:
  my $can_aggregate = $aggregator->disaggregate_types( $types_array_ref );
  # The idea is that if the aggregator knows how to create features of
  # a type in the list, it will add to the list any types that are
  # possible components of that type.

  # The return value of the disaggregate_types() method can be used to
  # determine whether or not the aggregate() method of this aggregator
  # can be used on features of the given types.
  my $should_aggregate = $aggregator->disaggregate_types( $types_array_ref );
  # Given $feature_provider, a CollectionProviderI of some sort,
  my $feature_collection =
    $feature_provider->get_collection( -types => $types_array_ref );
  if( $should_aggregate ) {
    $aggregator->aggregate( $feature_collection );
  }

=head1 DESCRIPTION

An AggregatorI is a factory for grouping features into components of a
containing (parent) feature.  Features of type 'exon' might be
aggregated into 'transcript's, for example.  The newly created
features B<contain> their component features, so if you use a
'transcript' aggregator to create a new feature of type 'transcript'
from a set of features of type 'exon', the new transcript feature will
contain the exon features (every feature is a
L<Bio::SeqFeature::CollectionI> of the component features that it
contains).

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

use Bio::Root::RootI;
@ISA = qw( Bio::Root::RootI );

use vars '$VERSION';
$VERSION = '1.00';

=head2 aggregate

 Title   : aggregate
 Usage   : my $did_aggregate =
             $aggregator->aggregate( $feature_array_ref [, $provider] )
 Function: Augment a list of L<Bio::SeqFeatureI> objects with any new
           containing (parent) features that can be created by this
           aggregator.
 Returns : A true value iff the given list was modified by the call.
 Args    : see below
 Status  : Public

  This method will create features that are aggregates of some or all
  of those features contained in the given list.  The new features, if
  there are any, will be added to the end of the list and the method
  will return true.  If no new features could be created the method
  will return a false value.

  The aggregate(..) method accepts its feature list by reference, and
  will modify the referent list directly.  It also optionally accepts
  a L<Bio::SeqFeature::CollectionProviderI> object; this is used by
  some aggregators to help determine the features that may be
  contained in a new parent.

=cut

sub aggregate {
  shift->throw_not_implemented();
}

=head2 disaggregate_types

 Title   : disaggregate_types
 Usage   : my $can_aggregate =
             $aggregator->disaggregate_types( $types_array_ref );
 Function: Augment a list of L<Bio::SeqFeature::TypeI> objects with
           possible components of those types.
 Returns : A true value iff this aggregator can aggregate any of the given
           types.
 Args    : A reference to a list of L<Bio::SeqFeature::TypeI> objects.
 Status  : Public

  The purpose of this method is both to determine whether this
  aggregator can be used to create features of the given types and
  also to ensure that the list contains all possible types that could
  be aggregated into those types.  Some examples: Let\'s say this
  aggregator can aggregate $foo_type and $bar_type features into
  $baz_type features, and $foo_type and $bozo_type features into
  $shezam_type features.  Then calling
  $aggregator->disaggregate_types( [ $baz_type, $bar_type ] ) will
  return true and cause the given list to become [ $baz_type,
  $bar_type, $foo_type ].  Calling $aggregator->disaggregate_types( [
  $bozo_type ] ) will return false and leave the list unchanged,
  because although the aggregator can use $bozo_type features to make
  $shezam_type features, the given list contained neither of the types
  that this aggregator can create.

=cut

sub disaggregate_types {
  shift->throw_not_implemented();
}

1;

__END__
