package Bio::DB::AggregatingFeatureProviderI;

# $Id $
# A Bio::DB::FeatureProviderI with additional methods to support the
# use of Bio::SeqFeature::AggregatorI objects for creating meta-features.

=head1 NAME

Bio::DB::AggregatingFeatureProviderI -- A provider of collections of
sequence features from a database or other non-trivial backing store,
with additional methods to support the use of
L<Bio::SeqFeature::AggregatorI>s.

=head1 SYNOPSIS

=head1 DESCRIPTION

The Bio::DB::AggregatingFeatureProviderI interface provides access to
Bio::SeqFeature::CollectionIs stored in a database or other (generally
external) backing store.  It is a Bio::SeqFeature::CollectionProviderI
with additional methods for managing the connection, caching, etc.,
and for maintaining a set of L<Bio::SeqFeature::AggregatorI>s to use
when creating features.

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
@ISA = qw( Bio::DB::FeatureProviderI );

use vars '$VERSION';
$VERSION = '1.00';

=head2 auto_aggregate

 Title   : auto_aggregate
 Usage   : my $auto_aggregate = $provider->auto_aggregate( [$new_value] );
 Function: Getter/Setter for the auto_aggregate flag.
 Returns : The current (or former, if used as a set method) boolean value
           of the auto_aggregate flag.
 Args    : [optional] a new value for the auto_aggregate flag.
 Status  : Public

   If auto_aggregate is true then feature collections returned by the
   get_collection(..) method will include any containing (parent)
   features that can be created using any aggregator in the
   aggregators() list, processed in the order that they appear in that
   list.  This method gets/sets the auto_aggregate flag.

=cut

sub auto_aggregate {
  shift->throw_not_implemented();
}

=head2 add_aggregator

 Title   : add_aggregator
 Usage   : $provider->add_aggregator( $aggregator );
 Function: Add an aggregator to the end of the list of aggregators.
 Returns : nothing
 Args    : An L<Bio::SeqFeature::AggregatorI> object.
 Status  : Public

   If auto_aggregate() is true, the aggregators in the aggregators()
   list will be processed, in order, on any
   L<Bio::SeqFeature::CollectionI> before it is returned by the
   get_collection() method.  The add_aggregator method adds an
   aggregator to the end of this list.  If you would like to add an
   aggregator elsewhere in the list, or otherwise control the list\'s
   contents in a more powerful way, use the aggregators() method as a
   setter.

=cut

sub add_aggregator {
  shift->throw_not_implemented();
}

=head2 aggregators

 Title   : aggregators
 Usage   : my @aggregators = $provider->aggregators( [@new_aggregators] );
 Function: Getter/Setter for the aggregators list.
 Returns : The current (or former, if used as a set method) list of
           L<Bio::SeqFeature::AggregatorI> objects.
 Args    : [optional] a new value for the aggregators list.
 Status  : Public

   If auto_aggregate is true then feature collections returned by the
   get_collection(..) method will include any containing (parent)
   features that can be created using any aggregator in the
   aggregators() list, processed in the order that they appear in that
   list.  This method returns the list and optionally changes the list.

=cut

sub aggregators {
  shift->throw_not_implemented();
}

=head2 clear_aggregators

 Title   : clear_aggregators
 Usage   : my @former_aggregators = $provider->clear_aggregators();
 Function: Clears the aggregators list.
 Returns : The former list of L<Bio::SeqFeature::AggregatorI> objects.
 Args    : none
 Status  : Public

   If auto_aggregate is true then feature collections returned by the
   get_collection(..) method will include any containing (parent)
   features that can be created using any aggregator in the
   aggregators() list, processed in the order that they appear in that
   list.  This method makes that list empty, and returns its former
   contents.

=cut

sub clear_aggregators {
  shift->throw_not_implemented();
}

1;

__END__
