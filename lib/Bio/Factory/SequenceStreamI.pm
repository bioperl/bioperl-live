#
# BioPerl module for Bio::Factory::SequenceStreamI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::SequenceStreamI - Interface describing the basics of a Sequence Stream.

=head1 SYNOPSIS

  # get a SequenceStreamI object somehow like with SeqIO
  use Bio::SeqIO;
  my $in = Bio::SeqIO->new(-file => '< fastafile');
  while( my $seq = $in->next_seq ) {
  }

=head1 DESCRIPTION

This interface is for describing objects which produces
Bio::PrimarySeqI objects or processes Bio::PrimarySeqI objects to a
data stream.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Factory::SequenceStreamI;
use strict;

use base qw(Bio::Root::RootI);

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = stream->next_seq
 Function: Reads the next sequence object from the stream and returns it.

           Certain driver modules may encounter entries in the stream that
           are either misformatted or that use syntax not yet understood
           by the driver. If such an incident is recoverable, e.g., by
           dismissing a feature of a feature table or some other non-mandatory
           part of an entry, the driver will issue a warning. In the case
           of a non-recoverable situation an exception will be thrown.
           Do not assume that you can resume parsing the same stream after
           catching the exception. Note that you can always turn recoverable
           errors into exceptions by calling $stream->verbose(2).
 Returns : a Bio::Seq sequence object
 Args    : none

See L<Bio::Root::RootI>

=cut

sub next_seq {
    shift->throw_not_implemented();
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object

=cut

sub write_seq {
    shift->throw_not_implemented();
}

=head2 sequence_factory

 Title   : sequence_factory
 Usage   : $seqio->sequence_factory($seqfactory)
 Function: Get the Bio::Factory::SequenceFactoryI
 Returns : Bio::Factory::SequenceFactoryI
 Args    : none


=cut

sub sequence_factory{
    shift->throw_not_implemented();
}

1;
