#
# BioPerl module for Bio::Seq::RichSeqI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::RichSeqI - interface for sequences from rich data sources, mostly databases

=head1 SYNOPSIS

    @secondary   = $richseq->get_secondary_accessions;
    $division    = $richseq->division;
    $mol         = $richseq->molecule;
    @dates       = $richseq->get_dates;
    $seq_version = $richseq->seq_version;
    $pid         = $richseq->pid;
    @keywords    = $richseq->get_keywords;

=head1 DESCRIPTION

This interface extends the L<Bio::SeqI> interface to give additional
functionality to sequences with richer data sources, in particular from database
sequences (EMBL, GenBank and Swissprot). For a general implementation, please
see the documentation for L<Bio::Seq::RichSeq>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::RichSeqI;
use strict;

use base qw(Bio::SeqI);


=head2 get_secondary_accessions

 Title   : get_secondary_accessions
 Usage   : 
 Function: Get the secondary accessions for a sequence.

           An implementation that allows modification of this array
           property should provide the methods add_secondary_accession
           and remove_secondary_accessions, with obvious purpose.

 Example :
 Returns : an array of strings
 Args    : none


=cut

sub get_secondary_accessions{
   my ($self,@args) = @_;

   $self->throw("hit get_secondary_accessions in interface definition - error");

}


=head2 division

 Title   : division
 Usage   :
 Function: Get (and set, depending on the implementation) the divison for
           a sequence.

           Examples from GenBank are PLN (plants), PRI (primates), etc.
 Example :
 Returns : a string
 Args    :


=cut

sub division{
   my ($self,@args) = @_;

   $self->throw("hit division in interface definition - error");

}


=head2 molecule

 Title   : molecule
 Usage   :
 Function: Get (and set, depending on the implementation) the molecule
           type for the sequence.

           This is not necessarily the same as Bio::PrimarySeqI::alphabet(),
           because it is databank-specific.
 Example :
 Returns : a string
 Args    :


=cut

sub molecule{
   my ($self,@args) = @_;

   $self->throw("hit molecule in interface definition - error");
}

=head2 pid

 Title   : pid
 Usage   :
 Function: Get (and set, depending on the implementation) the PID property
           for the sequence.
 Example :
 Returns : a string
 Args    :


=cut

sub pid {
   my ($self,@args) = @_;

   $self->throw("hit pid in interface definition - error");
}

=head2 get_dates

 Title   : get_dates
 Usage   :
 Function: Get (and set, depending on the implementation) the dates the
           databank entry specified for the sequence

           An implementation that allows modification of this array
           property should provide the methods add_date and
           remove_dates, with obvious purpose.

 Example :
 Returns : an array of strings
 Args    :


=cut

sub get_dates{
   my ($self,@args) = @_;

   $self->throw("hit get_dates in interface definition - error");

}


=head2 seq_version

 Title   : seq_version
 Usage   :
 Function: Get (and set, depending on the implementation) the version string
           of the sequence.
 Example :
 Returns : a string
 Args    :
 Note    : this differs from Bio::PrimarySeq version() in that this explicitly
           refers to the sequence record version one would find in a typical
           sequence file.  It is up to the implementation whether this is set
           separately or falls back to the more generic Bio::Seq::version()

=cut

sub seq_version{
   my ($self,@args) = @_;

   $self->throw("hit seq_version in interface definition - error");

}

=head2 get_keywords

 Title   : get_keywords
 Usage   : $obj->get_keywords()
 Function: Get the keywords for this sequence object.

           An implementation that allows modification of this array
           property should provide the methods add_keyword and
           remove_keywords, with obvious purpose.

 Returns : an array of strings
 Args    : 


=cut

sub get_keywords {
   my ($self) = @_;
   $self->throw("hit keywords in interface definition - error");
}

1;
