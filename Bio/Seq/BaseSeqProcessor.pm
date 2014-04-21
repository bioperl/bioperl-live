#
# BioPerl module for Bio::Seq::BaseSeqProcessor
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gmx.net, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
# 
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::BaseSeqProcessor - Base implementation for a SequenceProcessor

=head1 SYNOPSIS

    # you need to derive your own processor from this one

=head1 DESCRIPTION

This provides just a basic framework for implementations of
L<Bio::Factory::SequenceProcessorI>.

Essentially what it does is support a parameter to new() to set
sequence factory and source stream, and a next_seq() implementation
that will use a queue to be filled by a class overriding
process_seq().

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::BaseSeqProcessor;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::Root::Root Bio::Factory::SequenceProcessorI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Seq::BaseSeqProcessor->new();
 Function: Builds a new Bio::Seq::BaseSeqProcessor object 
 Returns : an instance of Bio::Seq::BaseSeqProcessor
 Args    : Named parameters. Currently supported are
             -seqfactory  the Bio::Factory::SequenceFactoryI object to use
             -source_stream the Bio::Factory::SequenceStreamI object to
                          which we are chained


=cut

sub new {
    my($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);

    my ($stream,$fact) =
	$self->_rearrange([qw(SOURCE_STREAM SEQFACTORY)], @args);

    $self->{'_queue'} = [];
    $self->sequence_factory($fact) if $fact;
    $self->source_stream($stream) if $stream;
    
    return $self;
}

=head1 L<Bio::Factory::SequenceProcessorI> methods

=cut

=head2 source_stream

 Title   : source_stream
 Usage   : $obj->source_stream($newval)
 Function: Get/set the source sequence stream for this sequence
           processor.

 Example : 
 Returns : A Bio::Factory::SequenceStreamI compliant object
 Args    : on set, new value (a Bio::Factory::SequenceStreamI compliant
           object)


=cut

sub source_stream{
    my $self = shift;

    if(@_) {
	my $stream = shift;
	my $fact = $stream->sequence_factory();
	$self->sequence_factory($fact)
	    unless $self->sequence_factory() || (! $fact);
	return $self->{'source_stream'} = $stream;
    }
    return $self->{'source_stream'};
}

=head1 L<Bio::Factory::SequenceStreamI> methods

=cut

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = stream->next_seq
 Function: Reads the next sequence object from the stream and returns it.

           This implementation will obtain objects from the source
           stream as necessary and pass them to process_seq() for
           processing. This method will return the objects one at a
           time that process_seq() returns.

 Returns : a Bio::Seq sequence object
 Args    : none

See L<Bio::Factory::SequenceStreamI::next_seq>

=cut

sub next_seq{
    my $self = shift;
    my $seq;

    # if the queue is empty, fetch next from source and process it
    if(@{$self->{'_queue'}} == 0) {
	my @seqs = ();
	while($seq = $self->source_stream->next_seq()) {
	    @seqs = $self->process_seq($seq);
	    # we may get zero seqs returned
	    last if @seqs;
	}
	push(@{$self->{'_queue'}}, @seqs) if @seqs;
    }
    # take next from the queue of seqs
    $seq = shift(@{$self->{'_queue'}});
    return $seq;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: Writes the result(s) of processing the sequence object into
           the stream.

           You need to override this method in order not to alter
           (process) sequence objects before output.

 Returns : 1 for success and 0 for error. The method stops attempting
           to write objects after the first error returned from the
           source stream. Otherwise the return value is the value
           returned from the source stream from writing the last
           object resulting from processing the last sequence object
           given as argument.

 Args    : Bio::SeqI object, or an array of such objects

=cut

sub write_seq{
    my ($self, @seqs) = @_;
    my $ret;
    foreach my $seq (@seqs) {
        foreach my $processed ($self->process_seq($seq)) {
            $ret = $self->source_stream->write_seq($seq);
            return unless $ret;
        }
    }
    return $ret;
}

=head2 sequence_factory

 Title   : sequence_factory
 Usage   : $seqio->sequence_factory($seqfactory)
 Function: Get the Bio::Factory::SequenceFactoryI
 Returns : Bio::Factory::SequenceFactoryI
 Args    : none


=cut

sub sequence_factory{
    my $self = shift;

    return $self->{'sequence_factory'} = shift if @_;
    return $self->{'sequence_factory'};
}

=head2 object_factory

 Title   : object_factory
 Usage   : $obj->object_factory($newval)
 Function: This is an alias to sequence_factory with a more generic name.
 Example : 
 Returns : a L<Bio::Factory::ObjectFactoryI> compliant object
 Args    : on set, new value (a L<Bio::Factory::ObjectFactoryI> 
           compliant object or undef, optional)


=cut

sub object_factory{
    return shift->sequence_factory(@_);
}

=head2 close

 Title   : close
 Usage   :
 Function: Closes the stream. We override this here in order to cascade
           to the source stream.
 Example :
 Returns : 
 Args    : none


=cut

sub close{
    my $self = shift;
    return $self->source_stream() ? $self->source_stream->close(@_) : 1;
}

=head1 To be overridden by a derived class

=cut

=head2 process_seq

 Title   : process_seq
 Usage   :
 Function: This is the method that is supposed to do the actual
           processing. It needs to be overridden to do what you want
           it to do.

           Generally, you do not have to override or implement any other
           method to derive your own sequence processor.

           The implementation provided here just returns the unaltered
           input sequence and hence is not very useful other than
           serving as a neutral default processor.

 Example :
 Returns : An array of zero or more Bio::PrimarySeqI (or derived
           interface) compliant object as the result of processing the
           input sequence.
 Args    : A Bio::PrimarySeqI (or derived interface) compliant object
           to be processed.


=cut

sub process_seq{
    my ($self,$seq) = @_;

    return ($seq);
}

1;
