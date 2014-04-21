#
# BioPerl module for Bio::Factory::SequenceProcessorI
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

Bio::Factory::SequenceProcessorI - Interface for chained sequence 
                                   processing algorithms

=head1 SYNOPSIS

    use Bio::SeqIO;
    use MySeqProcessor; # is-a Bio::Factory::SequenceProcessorI

    # obtain your source stream, e.g., an EMBL file
    my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'embl');
    # create your processor (it must implement this interface)
    my $seqalgo = MySeqProcessor->new();
    # chain together
    $seqalgo->source_stream($seqin);
    # you could create more processors and chain them one after another
    # ...
    # finally, the last link in the chain is your SeqIO stream
    my $seqpipe = $seqalgo;

    # once you've established the pipeline, proceed as if you had a
    # single SeqIO stream
    while(my $seq = $seqpipe->next_seq()) {
	# ... do something ...
    }

=head1 DESCRIPTION

This defines an interface that allows seamless chaining of sequence
processing algorithms encapsulated in modules while retaining the
overall Bio::SeqIO interface at the end of the pipeline.

This is especially useful if you want an easily configurable
processing pipeline of re-usable algorithms as building blocks instead
of (hard-)coding the whole algorithm in a single script.

There are literally no restrictions as to what an individual module
can do with a sequence object it obtains from the source stream before
it makes it available through its own next_seq() method. It can
manipulate the sequence object, but otherwise keep it intact, but it
can also create any number of new sequence objects from it, or it can
discard some, or any combination thereof. The only requirement is that
its next_seq() method return Bio::PrimarySeqI compliant objects. In
order to play nice, if a processor creates new objects it should try
to use the same sequence factory that the source stream uses, but this
is not strongly mandated.

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


package Bio::Factory::SequenceProcessorI;
use strict;
use Bio::Root::RootI;

use base qw(Bio::Factory::SequenceStreamI);

=head2 source_stream

 Title   : source_stream
 Usage   : $obj->source_stream($newval)
 Function: Get/set the source sequence stream for this sequence
           processor.

           An implementation is not required to allow set, but will
           usually do so.

 Example : 
 Returns : A Bio::Factory::SequenceStreamI compliant object
 Args    : on set, new value (a Bio::Factory::SequenceStreamI compliant
           object)


=cut

sub source_stream{
    shift->throw_not_implemented();
}

=head1 Bio::Factory::SequenceStreamI methods

 The requirement to implement these methods is inherited from
 L<Bio::Factory::SequenceStreamI>. An implementation may not
 necessarily have to implement all methods in a meaningful way. Which
 methods will be necessary very much depends on the context in which
 an implementation of this interface is used. E.g., if it is only used
 for post-processing sequences read from a SeqIO stream, write_seq()
 will not be used and hence does not need to be implemented in a
 meaningful way (it may in fact even throw an exception).

 Also, since an implementor will already receive built objects from a
 sequence stream, sequence_factory() may or may not be relevant,
 depending on whether the processing method does or does not involve
 creating new objects.

=cut

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = stream->next_seq
 Function: Reads the next sequence object from the stream and returns it.

           In the case of a non-recoverable situation an exception
           will be thrown.  Do not assume that you can resume parsing
           the same stream after catching the exception. Note that you
           can always turn recoverable errors into exceptions by
           calling $stream->verbose(2).

 Returns : a Bio::Seq sequence object
 Args    : none

See L<Bio::Root::RootI>

=cut

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object

=cut

=head2 sequence_factory

 Title   : sequence_factory
 Usage   : $seqio->sequence_factory($seqfactory)
 Function: Get the Bio::Factory::SequenceFactoryI
 Returns : Bio::Factory::SequenceFactoryI
 Args    : none


=cut

1;
