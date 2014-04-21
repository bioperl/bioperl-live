#
# BioPerl module for Bio::UpdateableSeqI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by David Block <dblock@gene.pbi.nrc.ca>
#
# Copyright David Block
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::UpdateableSeqI - Descendant of Bio::SeqI that allows updates

=head1 SYNOPSIS

See Bio::SeqI for most of the documentation.
See the documentation of the methods for further details.

=head1 DESCRIPTION

Bio::UpdateableSeqI is an interface for Sequence objects which are
expected to allow users to perform basic editing functions (update/delete)
on their component SeqFeatures.

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

=head1 AUTHOR - David Block

Email dblock@gene.pbi.nrc.ca

=head1 CONTRIBUTORS

Ewan Birney forced me to this...

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::UpdateableSeqI;
use strict;
use Carp;

# Object preamble - inherits from Bio::Root::Root



use base qw(Bio::SeqI);


=head2 delete_feature

 Title   : delete_feature
 Usage   : my $orphanlist=$self->delete_feature($feature,$transcript,$gene);
 Function: deletes the specified $feature from the given transcript, if $transcript is sent and exists and $feature is a feature of $transcript,
           or from $gene if the $feature is a feature of $gene, or from $self if $transcript and $gene are not sent.  Keeps track of the features
           of the $gene object that may be left as orphans and returns them as a listref.
 Example : I want to delete transcript 'abc' of gene 'def', with three exons, leaving only transcript 'ghi' with two exons.
           This will leave exons 1 and 3 part of 'ghi', but exon 2 will become an orphan.
           my $orphanlist=$seq->delete_feature($transcript{'abc'},undef,$gene{'def'});
           $orphanlist is a reference to a list containing $exon{'2'};
 Returns : a listref of orphaned features after the deletion of $feature (optional)
 Args    : $feature - the feature to be deleted
           $transcript - the transcript containing the $feature, so that a $feature can be removed from only one transcript when there are multiple
                         transcripts in a gene.
           $gene - the gene containing the $transcript and/or the $feature


=cut

sub delete_feature{
   my ($self,$feature,$transcript,$gene) = @_;

   $self->throw_not_implemented();
}


1;
