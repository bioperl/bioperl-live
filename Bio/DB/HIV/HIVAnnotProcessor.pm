# $Id: HIVAnnotProcessor.pm 221 2008-12-11 13:05:24Z maj $
#
# BioPerl module for HIVAnnotProcessor
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Mark A. Jensen <maj@fortinbras.us>
#
# Copyright Mark A. Jensen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

HIVAnnotProcessor - Adds HIV-specific annotations to Bio::SeqIO streams

=head1 SYNOPSIS

   sub get_Stream_by_query {
       my ($self, $query ) = @_;
       my $stream = $self->get_seq_stream('-query' => $query, '-mode'=>'query');
       return new Bio::DB::HIV::HIVAnnotProcessor( -hiv_query=>$query, 
                                                   -source_stream=>$stream );
   }
 

=head1 DESCRIPTION

Bio::DB::HIV::HIVAnnotProcessor is chained to the C<next_seq> of a sequence stream returned from a query to the Los Alamos HIV sequence database made using L<Bio::DB::HIV> and L<Bio::DB::Query::HIVQuery>. It adds the annotations obtained in the C<Bio::DB::Query::HIVQuery> to the Bio::Seq objects themselves via the C<$seq-E<gt>annotation> method.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mark A. Jensen

Email maj@fortinbras.us

=head1 CONTRIBUTORS

Mark A. Jensen

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::HIV::HIVAnnotProcessor;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use base qw( Bio::Root::Root);

=head1 Constructor

=head2 new

 Title   : new
 Usage   : my $obj = new HIVAnnotProcessor();
 Function: Builds a new HIVAnnotProcessor object 
 Returns : an instance of HIVAnnotProcessor
 Args    :

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($hiv_query, $source_stream) =
	$self->_rearrange([qw(HIV_QUERY SOURCE_STREAM)], @args); 


    $hiv_query                  && $self->hiv_query($hiv_query);
    $source_stream              && $self->source_stream($source_stream);

    return $self;
}

=head1 Bio::Factory::SequenceProcessorI compliance

=head2 source_stream

 Title   : source_stream
 Usage   : $hap->source_stream($newval)
 Function: 
 Example : 
 Returns : value of source_stream (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub source_stream{
    my $self = shift;
    if (@_) {
	$self->throw(-class=>'Bio::Root::BadParameter',
		     -text=>'Requires a Bio::SeqIO as argument',
		     -value=>$_[0]) unless $_[0]->isa('Bio::SeqIO');
    }
    return $self->{'source_stream'} = shift if @_;
    return $self->{'source_stream'};
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seqobj = stream->next_seq
 Function: Reads the next sequence object from the stream, 
         : adds annotations from the HIVQuery object according
         : to the sequence id, and returns sequence object
 Returns : a Bio::Seq sequence object
 Args    : none

=cut

sub next_seq {
    my $self = shift;
    my $q = $self->hiv_query;
    my $seqo = $self->source_stream->next_seq;
    return $seqo unless ($q && $seqo);
    
    my $ac = $q->get_annotations_by_id($seqo->primary_id);
    $seqo->annotation($ac) if $ac;
    my $acc = $q->get_accessions_by_id($seqo->primary_id);
    $seqo->accession_number($acc) if $acc;

    return $seqo;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $seqobj->write_seq
 Function: for HIVAnnotProcessor, throw an exception
 Example :
 Returns : Bio::Root::IOException
 Args    :

=cut

sub write_seq{
   my ($self,@args) = @_;
   $self->throw(-class=>'Bio::Root::IOException',
		-text=>'This stream is read-only',
		-value=>"");
}

=head1 HIVAnnotProcessor-specific methods

=head2 hiv_query

 Title   : hiv_query
 Usage   : $obj->hiv_query($newval)
 Function: 
 Example : 
 Returns : value of hiv_query (a Bio::DB::Query::HIVQuery object)
 Args    : on set, new value (an HIVQuery object, optional)

=cut

sub hiv_query{
    my $self = shift;
    if (@_) {
	$self->throw(-class=>'Bio::Root::BadParameter',
		     -text=>'Requires a Bio::DB::Query::HIVQuery as argument',
		     -value=>$_[0]) unless ref $_[0] && $_[0]->isa('Bio::DB::Query::HIVQuery');
    }
    return $self->{'hiv_query'} = shift if @_;
    return $self->{'hiv_query'};
}


1;

