
#
# BioPerl module for Bio::SeqI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqI - Abstract Interface of Sequence (with features)

=head1 SYNOPSIS

    # get a Bio::SeqI somehow, eg, via SeqIO. See SeqIO info
    # for more information

    $seqio  = Bio::SeqIO->new ( '-format' => 'Fasta' , -file => 'myfile.fasta');
    $seqobj = $seqio->next_seq();

    # features must implement Bio::SeqFeatureI

    @features = $seqobj->top_SeqFeatures(); # just top level
    @features = $seqobj->all_SeqFeatures(); # descend into sub features

    $seq      = $seqobj->seq(); # actual sequence as a string
    $seqstr   = $seqobj->subseq(10,50);    
    $ann      = $seqobj->annotation(); # annotation object

=head1 DESCRIPTION

SeqI is the abstract interface of annotated Sequence. These methods
are those which you can be guarenteed to get for any annseq. There aren't
many here, because too many complicated functions here prevent implementations
which are just wrappers around a database or similar delayed mechanisms.

Most of the clever stuff happens inside the SeqFeatureI system.

A good reference implementation is Bio::Seq which is a pure perl
implementation of this class with alot of extra pieces for extra
manipulation.  However, if you want to be able to use any sequence
object in your analysis, if you can do it just using these methods,
then you know you will be future proof and compatible with other
implementations of Seq.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::SeqI;
use strict;

use vars qw(@ISA);
use Bio::PrimarySeqI;
use Carp;

# Object preamble - inheriets from Bio::PrimarySeqI

@ISA = qw(Bio::PrimarySeqI);

=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   : my @features = $seq->top_SeqFeatures();
 Function: (may be deprecated in the future)
 Returns : An array of top level sequence features, 
           features which contain features are NOT flattened
 Args    : none


=cut

sub top_SeqFeatures{
   my ($self) = @_;
   $self->throw_not_implemented;
}


=head2 all_SeqFeatures

 Title   : all_SeqFeatures
 Usage   : @features = $annseq->all_SeqFeatures()
 Function: returns all SeqFeatures, included sub SeqFeatures
 Returns : an array
 Args    : none


=cut

sub all_SeqFeatures{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 seq

 Title   : seq
 Usage   : my $sequence  = $seq->seq();
 Function: Get/Set the sequence as a string
 Returns : string
 Args    : (optional) sequence string to set 
           (if permitted by the implementation)

=cut

sub seq{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 write_GFF

 Title   : write_GFF
 Usage   : $seq->write_GFF(\*OUTPUTHANDLE);
 Function: Utility methods which creates GFF output
 Returns : none
 Args    : (optional) filehandle to write to (default is STDOUT)


=cut

sub write_GFF{
   my ($self,$fh) = @_;

   $fh || do { $fh = \*STDOUT; };

   foreach my $sf ( $self->all_SeqFeatures() ) {
       print $fh $sf->gff_string, "\n";
   }
}

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($seq_obj)
 Function: Get/Set Annotation object associated with the sequence.  
           This is where features with no explicit location are aggregated. 
 Returns : Bio::AnnotationCollectionI object
 Args    : None or Bio::Annotation object

See L<Bio::AnnotationCollectionI> and L<Bio::Annotation::Collection>
for more information

=cut

sub annotation {
   my ($obj) = @_;
   $obj->throw_not_implemented();
}

=head2 feature_count

 Title   : feature_count
 Usage   : $seq->feature_count()
 Function: Return the number of SeqFeatures attached to a sequence
 Returns : number of SeqFeatures
 Args    : none


=cut

sub feature_count {
    my ($obj) = @_;
    $obj->throw_not_implemented();
}

=head2 species

 Title   : species
 Usage   : $species = $self->species();
 Function: Gets or sets the species
 Returns : Bio::Species object
 Args    : Bio::Species object or none;


=cut

sub species {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 primary_seq

 Title   : primary_seq
 Usage   : $obj->primary_seq($newval)
 Function: Get/Set the underlying primary sequence for a SeqI.
           This cannot always be assumed some implementations
           may not provide a primary sequence underneath.
 Returns : value of primary_seq
 Args    : newvalue (optional)

=cut

sub primary_seq {
    my ($self) = @_;
    $self->throw_not_implemented;
}

1;
