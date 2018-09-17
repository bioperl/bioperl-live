#
# BioPerl module for Bio::SeqFeature::Gene::TranscriptI
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

Bio::SeqFeature::Gene::TranscriptI - Interface for a feature representing a
         transcript of exons, promoter(s), UTR, and a poly-adenylation site.

=head1 SYNOPSIS

  #documentation needed

=head1 DESCRIPTION

A feature representing a transcript.

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

Email hlapp@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene::TranscriptI;
use strict;

use Carp;

use base qw(Bio::SeqFeatureI);

=head2 promoters

 Title   : promoters()
 Usage   : @proms = $transcript->promoters();
 Function: Get the promoter features of this transcript.

           Note that OO-modeling of regulatory elements is not stable yet.
           This means that this method might change or even disappear in a
           future release. Be aware of this if you use it.

 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           promoter regions or sites.
 Args    : 

=cut

sub promoters {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 exons

 Title   : exons()
 Usage   : @exons = $transcript->exons();
           @inital = $transcript->exons('Initial');
 Function: Get the individual exons this transcript comprises of, or all exons
           of a specified type.

           Refer to the documentation of the class that produced this 
           transcript object for information about the possible types.

           See Bio::SeqFeature::Gene::ExonI for properties of the
           returned objects.

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects
 Args    : An optional string specifying the type of the exon.

=cut

sub exons {
    my ($self, $type) = @_;
    $self->throw_not_implemented();
}

=head2 introns

 Title   : introns()
 Usage   : @introns = $transcript->introns();
 Function: Get all introns this transcript comprises of.
 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           introns.
 Args    : 


=cut

sub introns {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 poly_A_site

 Title   : poly_A_site()
 Usage   : $polyAsite = $transcript->poly_A_site();
 Function: Get the poly-adenylation site of this transcript.
 Returns : A Bio::SeqFeatureI implementing object.
 Args    : 


=cut

sub poly_A_site {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 utrs

 Title   : utrs()
 Usage   : @utr_sites = $transcript->utrs();
 Function: Get the UTR regions this transcript comprises of.

           See Bio::SeqFeature::Gene::ExonI for properties of the
           returned objects.

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects
 Args    : 


=cut

sub utrs {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 mrna

 Title   : mrna()
 Usage   : $mrna = $transcript->mrna();
 Function: Get the mRNA of the transcript as a sequence object.

 Returns : A Bio::PrimarySeqI implementing object.
 Args    : 


=cut

sub mrna {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 cds

 Title   : cds()
 Usage   : $cds = $transcript->cds();
 Function: Get the CDS (coding sequence) of the transcript as a sequence
           object.

 Returns : A Bio::PrimarySeqI implementing object.
 Args    : 


=cut

sub cds {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 protein

 Title   : protein()
 Usage   : $protein = $transcript->protein();
 Function: Get the protein encoded by the transcript as a sequence object.

 Returns : A Bio::PrimarySeqI implementing object.
 Args    : 


=cut

sub protein {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 parent

 Title   : parent
 Usage   : $obj->parent($newval)
 Function: get the parent gene of the transcript
 Returns : value of parent - a Bio::SeqFeature::Gene::GeneStructureI-compliant object
 Args    : a Bio::SeqFeature::Gene::GeneStructureI-compliant object (optional)


=cut

sub parent{
   my $self = shift;
   if( @_ ) {
       my $value = shift;
       # I really mean ! defined $value - 
       # we will allow re-setting the parent to undef 
       if (! defined $value || 
	   $value->isa("Bio::SeqFeature::Gene::GeneStructureI")) {
	   $self->{'_parent'} = $value;
       } else {
	   $self->throw("$value must be a Bio::SeqFeature::Gene::GeneStructureI")
       }
    }
    return $self->{'_parent'};
}


1;






