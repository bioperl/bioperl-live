#
# BioPerl module for Bio::SeqI
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

Bio::SeqI - [Developers] Abstract Interface of Sequence (with features)

=head1 SYNOPSIS

    # Bio::SeqI is the interface class for sequences.

    # If you are a newcomer to bioperl, you should
    # start with Bio::Seq documentation. This
    # documentation is mainly for developers using
    # Bioperl.

    # Bio::SeqI implements Bio::PrimarySeqI
    $seq      = $seqobj->seq(); # actual sequence as a string
    $seqstr   = $seqobj->subseq(10,50);

    # Bio::SeqI has annotationcollections

    $ann      = $seqobj->annotation(); # annotation object

    # Bio::SeqI has sequence features
    # features must implement Bio::SeqFeatureI

    @features = $seqobj->get_SeqFeatures();     # just top level
    @features = $seqobj->get_all_SeqFeatures(); # descend into sub features

=head1 DESCRIPTION

Bio::SeqI is the abstract interface of annotated Sequences. These
methods are those which you can be guaranteed to get for any Bio::SeqI.
For most users of the package the documentation (and methods) in this
class are not at useful - this is a developers only class which
defines what methods have to be implemented by other Perl objects to
comply to the Bio::SeqI interface. Go "perldoc Bio::Seq" or "man
Bio::Seq" for more information.

There aren't many method here, because too many complicated functions here
would prevent implementations which are just wrappers around a database or
similar delayed mechanisms.

Most of the clever stuff happens inside the SeqFeatureI system.

A good reference implementation is Bio::Seq which is a pure perl
implementation of this class with a lot of extra pieces for extra
manipulation.  However, if you want to be able to use any sequence
object in your analysis, if you can do it just using these methods,
then you know you will be future proof and compatible with other
implementations of Seq.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::SeqI;
use strict;


# Object preamble - inherits from Bio::PrimarySeqI

use base qw(Bio::PrimarySeqI Bio::AnnotatableI Bio::FeatureHolderI);

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : my @feats = $seq->get_SeqFeatures();
 Function: retrieve just the toplevel sequence features attached to this seq
 Returns : array of Bio::SeqFeatureI objects
 Args    : none

This method comes through extension of Bio::FeatureHolderI. See
L<Bio::FeatureHolderI> and L<Bio::SeqFeatureI> for more information.

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : my @feats = $seq->get_all_SeqFeatures();
 Function: returns all SeqFeatures, including sub SeqFeatures
 Returns : an array of Bio::SeqFeatureI objects
 Args    : none

This method comes through extension of Bio::FeatureHolderI. See
L<Bio::FeatureHolderI> and L<Bio::SeqFeatureI> for more information.

=head2 feature_count

 Title   : feature_count
 Usage   : my $count = $seq->feature_count();
 Function: Return the number of SeqFeatures attached to a sequence
 Returns : integer representing the number of SeqFeatures
 Args    : none

This method comes through extension of Bio::FeatureHolderI. See
L<Bio::FeatureHolderI> for more information.

=head2 seq

 Title   : seq
 Usage   : my $string = $seq->seq();
 Function: Retrieves the sequence string for the sequence object
 Returns : string
 Args    : none


=cut

sub seq {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 write_GFF

 Title   : write_GFF
 Usage   : $seq->write_GFF(\*FILEHANDLE);
 Function: Convenience method to write out all the sequence features
           in GFF format to the provided filehandle (STDOUT by default)
 Returns : none
 Args    : [optional] filehandle to write to (default is STDOUT)


=cut

sub write_GFF {
   my ($self,$fh) = @_;

   $fh || do { $fh = \*STDOUT; };

   foreach my $sf ( $self->get_all_SeqFeatures() ) {
       print $fh $sf->gff_string, "\n";
   }

}

=head2 annotation

 Title   : annotation
 Usage   : my $ann = $seq->annotation($seq_obj);
 Function: retrieve the attached annotation object
 Returns : Bio::AnnotationCollectionI or none;

See L<Bio::AnnotationCollectionI> and L<Bio::Annotation::Collection>
for more information. This method comes through extension from
L<Bio::AnnotatableI>.

=head2 species

 Title   : species
 Usage   :
 Function: Gets or sets the species
 Example : my $species = $seq->species();
 Returns : Bio::Species object
 Args    : Bio::Species object or none;

See L<Bio::Species> for more information

=cut

sub species {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 primary_seq

 Title   : primary_seq
 Usage   : my $primaryseq = $seq->primary_seq($newval)
 Function: Retrieve the underlying Bio::PrimarySeqI object if available.
           This is in the event one has a sequence with lots of features
           but want to be able to narrow the object to just one with
           the basics of a sequence (no features or annotations).
 Returns : Bio::PrimarySeqI
 Args    : Bio::PrimarySeqI or none;

See L<Bio::PrimarySeqI> for more information

=cut

sub primary_seq {
    my ($self) = @_;
    $self->throw_not_implemented;
}

1;
