package Bio::SeqFeature::CompoundIterator;

# $Id$
# An iterator over Bio::SeqFeatureI objects that can contain other iterators.

=head1 NAME

Bio::SeqFeature::CompoundIterator - An iterator over Bio::SeqFeatureI objects.

=head1 SYNOPSIS

  my $iterator = $seq_feature_collection->iterator( @args );
  while( $iterator->has_more_elements() ) {
    my $seq_feature = $iterator->next_feature();
    # do something with the features...
  }

=head1 DESCRIPTION

  A compound iterator over Bio::SeqFeatureI objects.  This is an
  autodeleting list, aka a stream.  Reading from the iterator (using
  $iterator->next_feature) removes the returned feature from the
  iterator.  An iterator may return new objects or it may return
  preexisting objects.  In some circumstances it may return the same
  object twice.  See Bio::SeqFeature::CollectionI.  The compound
  iterator can iterate over features from other iterators after
  iterating over its own.

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
use Bio::SeqFeature::SimpleIterator;
use vars qw( $VERSION @ISA );

$VERSION = '0.01';
@ISA = qw( Bio::SeqFeature::SimpleIterator );

=head2 new

 Title   : new
 Usage   : $iterator = new Bio::SeqFeature::CompoundIterator( $feature_list_ref )
           OR
           $iterator = new Bio::SeqFeature::CompoundIterator( @feature_list )
 Function: Instantiates a new iterator over the given features
 Returns : a new Bio::SeqFeature::CompoundIterator
 Args    : A list (or list ref) of SeqFeatureI objects
 Status  : Public

=cut

sub new {
  my $caller = shift;

  my $self = $caller->SUPER::new( @_ );

  $self->{ '_next_iterators' } = [];

  return $self;
} # new(..)

## TODO: Document this.
sub add_next_iterator {
  my $self = shift;

  foreach my $next_iterator ( @_ ) {
    push( @{ $self->{ '_next_iterators' } }, $next_iterator );
  }
  return;
} # add_next_iterator(..)

## TODO: Document this.
sub get_next_iterators {
  my $self = shift;

  return @{ $self->{ '_next_iterators' } };
} # get_next_iterators()

=head2 next_feature

 Title   : next_feature
 Usage   : $seq_feature = $iterator->next_feature()
 Function: returns and removes the next feature from this iterator
 Returns : a Bio::SeqFeatureI, or undef if there are no more
 Args    : none
 Status  : Public

=cut

sub next_feature {
  my $self = shift;

  my $current_iterator = $self->{ '_current_iterator' };
  unless( $current_iterator ) {
    if( $self->{ '_started' } ) {
      return undef; # No more.
    } else {
      $self->{ '_started' } = 1;
      $current_iterator = $self->{ '_current_iterator' } = $self;
    }
  }
  my $next;
  while( defined( $current_iterator ) ) {
    if( $current_iterator == $self ) {
      if( $self->SUPER::has_more_features() ) {
        $next = $self->SUPER::next_feature();
      }
    } elsif( $current_iterator->has_more_features() ) {
      $next = $current_iterator->next_feature();
    }
    if( defined( $next ) ) {
      return $next;
    }
    $current_iterator =
      $self->{ '_current_iterator' } =
      shift( @{ $self->{ '_next_iterators' } } );
  }
  return undef; # No more.
} # next_feature()

=head2 has_more_features

 Title   : has_more_features
 Usage   : while( $iterator->has_more_features() ) { do something }
 Function: returns true iff there are features in this iterator
 Returns : true iff there are more features in this iterator
 Args    : none
 Status  : Public

=cut

sub has_more_features {
  my $self = shift;

  my $current_iterator = $self->{ '_current_iterator' };
  unless( $current_iterator ) {
    if( $self->{ '_started' } ) {
      return 0; # No more.
    } else {
      $self->{ '_started' } = 1;
      $current_iterator = $self->{ '_current_iterator' } = $self;
    }
  }
  my $has_next;
  while( defined( $current_iterator ) ) {
    if( $current_iterator == $self ) {
      $has_next = $self->SUPER::has_more_features();
    } else {
      $has_next = $current_iterator->has_more_features();
    }
    if( $has_next ) {
      return $has_next;
    }
    $current_iterator =
      $self->{ '_current_iterator' } =
      shift( @{ $self->{ '_next_iterators' } } );
  }
  return 0; # No more.
} # has_more_features()

1;

__END__
