# $Id$
#
# BioPerl module for Bio::SeqFeature::Gene::GeneStructureI
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Gene::GeneStructure - A feature representing an arbitrarily
           complex structure of a gene

=head1 SYNOPSIS

  #documentaion needed

=head1 DESCRIPTION

A feature representing a gene structure.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene::GeneStructureI;
use vars qw(@ISA);
use strict;

use Carp;
use Bio::SeqFeatureI;

@ISA = qw(Bio::SeqFeatureI);

=head2 transcripts

 Title   : transcripts()
 Usage   : @transcripts = $gene->transcripts();
 Function: Get the transcript features/sites of this gene structure.

           See Bio::SeqFeature::Gene::TranscriptI for properties of the
           returned objects.

 Returns : An array of Bio::SeqFeature::Gene::TranscriptI implementing objects
           representing the promoter regions or sites.
 Args    : 


=cut

sub transcripts {
    my ($self) = @_;

    $self->_abstractDeath('transcripts');
}

=head2 promoters

 Title   : promoters()
 Usage   : @prom_sites = $gene->promoters();
 Function: Get the promoter features/sites of this gene structure.

           Note that OO-modeling of regulatory elements is not stable yet.
           This means that this method might change or even disappear in a
           future release. Be aware of this if you use it.

 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           promoter regions or sites.
 Args    : 

=cut

sub promoters {
    my ($self) = @_;

    $self->_abstractDeath('promotors');
}

=head2 exons

 Title   : exons()
 Usage   : @exons = $gene->exons();
           @inital = $gene->exons('Initial');
 Function: Get all exon features or all exons of specified type of this gene
           structure.

           Refer to the documentation of the class that produced this gene
           structure object for information about the possible types.

           See Bio::SeqFeature::Gene::ExonI for properties of the
           returned objects.

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects
           representing the exon regions.
 Args    : An optional string specifying the type of the exon.

=cut

sub exons {
    my ($self, $type) = @_;

    $self->_abstractDeath('exons');
}

=head2 introns

 Title   : introns()
 Usage   : @introns = $gene->introns();
 Function: Get all introns of this gene structure.
 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           introns.
 Args    : 


=cut

sub introns {
    my ($self) = @_;

    $self->_abstractDeath('introns');    
}

=head2 poly_A_sites

 Title   : poly_A_sites()
 Usage   : @polyAsites = $gene->poly_A_sites();
 Function: Get the poly-adenylation features/sites of this gene structure.
 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           poly-adenylation regions or sites.
 Args    : 


=cut

sub poly_A_sites {
    my ($self) = @_;

    $self->_abstractDeath('poly_A_sites');
}

=head2 utrs

 Title   : utrs()
 Usage   : @utr_sites = $gene->utrs();
 Function: Get the UTR features/sites of this gene structure.

           See Bio::SeqFeature::Gene::ExonI for properties of the
           returned objects.

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects
           representing the UTR regions or sites.
 Args    : 


=cut

sub utrs {
    my ($self) = @_;

    $self->_abstractDeath('utrs');
}

1;




