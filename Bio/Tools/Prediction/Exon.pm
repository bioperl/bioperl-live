#
# BioPerl module for Bio::Tools::Prediction::Exon
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

Bio::Tools::Prediction::Exon - A predicted exon feature

=head1 SYNOPSIS

  # See documentation of methods.

=head1 DESCRIPTION

A feature representing a predicted exon. This class actually inherits
off Bio::SeqFeature::Gene::Exon and therefore has all that
functionality (also implements Bio::SeqFeatureI), plus a few methods
supporting predicted features, like various scores and a
significance. Even though these were inspired by GenScan results, at
least a subset should be generally useable for exon prediction
results.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp-at-gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Prediction::Exon;
use strict;


use base qw(Bio::SeqFeature::Gene::Exon);

sub new {
    my($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);

    return $self;
}


=head2 predicted_cds

 Title   : predicted_cds
 Usage   : $predicted_cds_dna = $exon->predicted_cds();
           $exon->predicted_cds($predicted_cds_dna);
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
 Usage   : $predicted_protein_seq = $exon->predicted_protein();
           $exon->predicted_protein($predicted_protein_seq);
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

=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function: 
 Returns : 
 Args    : 


=cut

sub significance {
    return shift->_tag_value('signif', @_);
}

=head2 start_signal_score

 Title   : start_signal_score
 Usage   : $sc = $obj->start_signal_score();
           $obj->start_signal_score($evalue);
 Function: Get/Set a score for the exon start signal (acceptor splice site
           or initiation signal).
 Returns : 
 Args    : 


=cut

sub start_signal_score {
    return shift->_tag_value('AccScore', @_);
}

=head2 end_signal_score

 Title   : end_signal_score
 Usage   : $sc = $obj->end_signal_score();
           $obj->end_signal_score($evalue);
 Function: Get/Set a score for the exon end signal (donor splice site
           or termination signal).
 Returns : 
 Args    : 


=cut

sub end_signal_score {
    return shift->_tag_value('DonScore', @_);
}

=head2 coding_signal_score

 Title   : coding_signal_score
 Usage   : $sc = $obj->coding_signal_score();
           $obj->coding_signal_score($evalue);
 Function: Get/Set a score for the exon coding signal (e.g., coding potential).
 Returns : 
 Args    : 


=cut

sub coding_signal_score {
    return shift->_tag_value('CodScore', @_);
}

#
# Everything else is just inherited from SeqFeature::Generic.
#

1;
