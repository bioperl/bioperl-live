# $Id $

=head1 NAME

Bio::SeqFeature::IteratorI - An iterator over Bio::SeqFeatureI objects.

=head1 SYNOPSIS

  my $iterator = $seq_feature_collection->iterator( @args );
  while( $iterator->has_more_elements() ) {
    my $seq_feature = $iterator->next_feature();
    # do something with the features...
  }

=head1 DESCRIPTION

  An iterator over Bio::SeqFeatureI objects.  This is an autodeleting
  list, aka a stream.  Reading from the iterator (using
  $iterator->next_feature) removes the returned feature from the
  iterator.  An iterator may return new objects or it may return
  preexisting objects.  In some circumstances it may return the same
  object twice (eg. if it is iterating over a list that contains that
  object twice).  See Bio::SeqFeature::CollectionI.

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

package Bio::SeqFeature::IteratorI;

use strict;
use vars qw( $VERSION @ISA );

$VERSION = '0.01';
@ISA = qw( Bio::Root::RootI );

=head2 next_feature

 Title   : next_feature
 Usage   : $seq_feature = $iterator->next_feature()
 Function: returns and removes the next feature from this iterator
 Returns : a Bio::SeqFeatureI, or undef if there are no more
 Args    : none
 Status  : Public

=cut

sub next_feature {
  shift->throw_not_implemented();
} # next_feature()

=head2 next_seq

 Title   : next_seq
 Usage   : $seq_feature = $iterator->next_seq()
 Function: returns and removes the next feature from this iterator
 Returns : a Bio::SeqFeatureI, or undef if there are no more
 Args    : none
 Status  : Public

  This method is a (glob ref) alias for next_feature.

=cut

  *next_seq = \&next_feature;

=head2 has_more_features

 Title   : has_more_features
 Usage   : while( $iterator->has_more_features() ) { do something }
 Function: returns true iff there are features in this iterator
 Returns : true iff there are more features in this iterator
 Args    : none
 Status  : Public

=cut

sub has_more_features {
  shift->throw_not_implemented();
} # has_more_features()

1;
