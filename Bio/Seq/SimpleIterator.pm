package Bio::Seq::SimpleIterator;

# $Id $
# An iterator over Bio::PrimarySeqI objects.

=head1 NAME

Bio::Seq::SimpleIterator - An iterator over Bio::PrimarySeqI objects.

=head1 SYNOPSIS

  my $iterator = $sequence_provider->get_Seq_stream( @args );
  while( $iterator->has_more_elements() ) {
    my $sequence = $iterator->next_seq();
    # do something with the sequences...
  }

=head1 DESCRIPTION

  An iterator over Bio::PrimarySeqI objects.  This is an autodeleting
  list, aka a stream.  Reading from the iterator (using
  $iterator->next_sequence) removes the returned sequence from the
  iterator.  An iterator may return new objects or it may return
  preexisting objects.  In some circumstances it may return the same
  object twice.  See L<Bio::DB::SequenceProviderI>.

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
use Bio::Root::Root;
use Bio::Seq::IteratorI;
use vars qw( $VERSION @ISA );

$VERSION = '0.01';
@ISA = qw( Bio::Root::Root Bio::Seq::IteratorI );

=head2 new

 Title   : new
 Usage   : $iterator = new Bio::Seq::SimpleIterator( $sequence_list_ref )
           OR
           $iterator = new Bio::Seq::SimpleIterator( @sequence_list )
 Function: Instantiates a new iterator over the given sequences
 Returns : a new Bio::Seq::SimpleIterator
 Args    : A list (or list ref) of L<Bio::PrimarySeqI> objects
 Status  : Public

=cut

sub new {
  my $class = shift;
  $class = ref($class) if ref($class);

  my $sequences;
  if( scalar( @_ ) > 1 ) {
    $sequences = [];
    push( @$sequences, @_ );
  } elsif( scalar( @_ ) == 1 ) {
    if( ref $_[ 0 ] eq 'ARRAY' ) {
      $sequences = shift;
    } else {
      $sequences = [ shift ];
    }
  }
  return bless {
                 '_sequences'  => $sequences
	       },$class;
} # new(..)

=head2 next_seq

 Title   : next_seq
 Usage   : $sequence = $iterator->next_seq()
 Function: returns and removes the next sequence from this iterator
 Returns : a L<Bio::PrimarySeqI>, or undef if there are no more
 Args    : none
 Status  : Public

=cut

sub next_seq {
  my $self = shift;

  if( $self->{ '_sequences' } ) {
    return shift @{ $self->{ '_sequences' } };
  }
  return undef;
} # next_sequence()

=head2 has_more_sequences

 Title   : has_more_sequences
 Usage   : while( $iterator->has_more_sequences() ) { do something }
 Function: returns true iff there are sequences in this iterator
 Returns : true iff there are more sequences in this iterator
 Args    : none
 Status  : Public

=cut

sub has_more_sequences {
  my $self = shift;
  return ( $self->{ '_sequences' } && scalar( @{ $self->{ '_sequences' } } ) );
} # has_more_sequences()

1;

__END__
