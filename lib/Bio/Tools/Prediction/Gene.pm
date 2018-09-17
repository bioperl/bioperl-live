#
# BioPerl module for Bio::Tools::Prediction::Gene
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Prediction::Gene - a predicted gene structure feature

=head1 SYNOPSIS

  #See documentation of methods.

=head1 DESCRIPTION

A feature representing a predicted gene structure. This class actually
inherits off Bio::SeqFeature::Gene::Transcript and therefore has all that
functionality, plus a few methods supporting predicted sequence features,
like a predicted CDS and a predicted translation.

Exons held by an instance of this class will usually be instances of
Bio::Tools::Prediction::Exon, although they do not have to be. Refer to the
documentation of the class that produced the instance.

Normally, you will not want to create an instance of this class yourself.
Instead, classes representing the results of gene structure prediction
programs will do that.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp-at-gmx.net or hilmar.lapp-at-pharma.novartis.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Prediction::Gene;
use strict;



use base qw(Bio::SeqFeature::Gene::Transcript);

sub new {
    my($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);

    my ($primary,$ptag) = $self->_rearrange([qw(PRIMARY PRIMARY_TAG)],@args);
    $self->primary_tag('predicted_gene') unless $primary || $ptag;

    return $self; 
}


=head2 predicted_cds

 Title   : predicted_cds
 Usage   : $predicted_cds_dna = $gene->predicted_cds();
           $gene->predicted_cds($predicted_cds_dna);
 Function: Get/Set the CDS (coding sequence) as predicted by a program.

           This method is independent of an attached_seq. There is no
           guarantee whatsoever that the returned CDS has anything to do
           (e.g., matches) with the sequence covered by the exons as annotated
           through this object.

 Example :
 Returns : A Bio::PrimarySeqI implementing object holding the DNA sequence
           defined as coding by a prediction of a program.
 Args    : On set, a Bio::PrimarySeqI implementing object holding the DNA 
           sequence defined as coding by a prediction of a program.

=cut

sub predicted_cds {
    my ($self, $cds) = @_;

    if(defined($cds)) {
	$self->{'_predicted_cds'} = $cds;
    }
    return $self->{'_predicted_cds'};
}

=head2 predicted_protein

 Title   : predicted_protein
 Usage   : $predicted_protein_seq = $gene->predicted_protein();
           $gene->predicted_protein($predicted_protein_seq);
 Function: Get/Set the protein translation as predicted by a program.

           This method is independent of an attached_seq. There is no
           guarantee whatsoever that the returned translation has anything to
           do with the sequence covered by the exons as annotated
           through this object, or the sequence returned by predicted_cds(),
           although it should usually be just the standard translation.

 Example :
 Returns : A Bio::PrimarySeqI implementing object holding the protein 
           translation as predicted by a program.
 Args    : On set, a Bio::PrimarySeqI implementing object holding the protein 
           translation as predicted by a program.

=cut

sub predicted_protein {
    my ($self, $aa) = @_;

    if(defined($aa)) {
	$self->{'_predicted_aa'} = $aa;
    }
    return $self->{'_predicted_aa'};
}

#
# Everything else is just inherited from SeqFeature::GeneStructure.
#

1;
