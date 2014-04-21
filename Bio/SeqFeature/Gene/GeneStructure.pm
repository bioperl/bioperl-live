#
# BioPerl module for Bio::SeqFeature::Gene::GeneStructure
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

Bio::SeqFeature::Gene::GeneStructure - A feature representing an arbitrarily complex structure of a gene

=head1 SYNOPSIS

  # See documentation of methods.

=head1 DESCRIPTION

A feature representing a gene structure. As of now, a gene structure
really is only a collection of transcripts. See
L<Bio::SeqFeature::Gene::TranscriptI> (interface) and
L<Bio::SeqFeature::Gene::Transcript> (implementation) for the features
of such objects.

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
 the bugs and their resolution. Bug reports can be submitted via the
 web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp-at-gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene::GeneStructure;
use vars qw($WeakRefs);
use strict;

BEGIN {
    eval "use Scalar::Util qw(weaken);";
    if ($@) {
	$Bio::SeqFeature::Gene::GeneStructure::WeakRefs = 0;  
    } else { $Bio::SeqFeature::Gene::GeneStructure::WeakRefs = 1; }
}


use base qw(Bio::SeqFeature::Generic Bio::SeqFeature::Gene::GeneStructureI);


sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);
    $self->_register_for_cleanup(\&gene_cleanup);
    my ($primary) =
	$self->_rearrange([qw(PRIMARY
			      )],@args);
    
    $primary = 'genestructure' unless $primary;
    $self->primary_tag($primary);
    $self->strand(0) if(! defined($self->strand()));
    return $self;
}

=head2 transcripts

 Title   : transcripts
 Usage   : @transcripts = $gene->transcripts();
 Function: Get the transcripts of this gene structure. Many gene structures
           will have only one transcript.

 Returns : An array of Bio::SeqFeature::Gene::TranscriptI implementing objects.
 Args    : 


=cut

sub transcripts {
    return @{shift->{'_transcripts'} || []};
}

=head2 add_transcript

 Title   : add_transcript()
 Usage   : $gene->add_transcript($transcript);
 Function: Add a transcript to this gene structure.
 Returns : 
 Args    : A Bio::SeqFeature::Gene::TranscriptI implementing object.


=cut

sub add_transcript {
    my ($self, $fea) = @_;

    if(!$fea || ! $fea->isa('Bio::SeqFeature::Gene::TranscriptI') ) {
	$self->throw("$fea does not implement Bio::SeqFeature::Gene::TranscriptI");
    }
    unless( exists $self->{'_transcripts'}  ) {
	$self->{'_transcripts'} = [];
    }
    $self->_expand_region($fea);
    if( $Bio::SeqFeature::Gene::GeneStructure::WeakRefs ) {
	$fea->parent(weaken $self);
    } else {
	$fea->parent($self);
    }
    push(@{$self->{'_transcripts'}}, $fea);
}

=head2 flush_transcripts

 Title   : flush_transcripts()
 Usage   : $gene->flush_transcripts();
 Function: Remove all transcripts from this gene structure.
 Returns : 
 Args    :


=cut

sub flush_transcripts {
    my ($self) = @_;    
    if( defined $self->{'_transcripts'} ) {
	foreach my $t ( grep {defined} @{$self->{'_transcripts'} || []} ) {
	    $t->parent(undef); # remove bkwds pointers
	    $t = undef;
	}
	delete($self->{'_transcripts'});	
    }
}

=head2 add_transcript_as_features

 Title   : add_transcript_as_features
 Usage   : $gene->add_transcript_as_features(@featurelist);
 Function: take a list of Bio::SeqFeatureI objects and turn them into a
           Bio::SeqFeature::Gene::Transcript object.  Add that transcript to the gene.
 Returns : nothing
 Args    : a list of Bio::SeqFeatureI compliant objects


=cut

sub add_transcript_as_features {
    my ($self,@features) = @_;
    my $transcript=Bio::SeqFeature::Gene::Transcript->new;
    foreach my $fea (@features) {
	if ($fea->primary_tag =~ /utr/i) { #UTR / utr/ 3' utr / utr5 etc.
	    $transcript->add_utr($fea);
	} elsif ($fea->primary_tag =~ /promot/i) { #allow for spelling differences
	    $transcript->add_promoter($fea);
	} elsif ($fea->primary_tag =~ /poly.*A/i) { #polyA, POLY_A, etc.
	    $transcript->poly_A_site($fea);
	} else {		#assume the rest are exons
	    $transcript->add_exon($fea);
	}
    }
    $self->add_transcript($transcript);
}


=head2 promoters

 Title   : promoters
 Usage   : @prom_sites = $gene->promoters();
 Function: Get the promoter features of this gene structure.

           This method basically merges the promoters returned by transcripts.

           Note that OO-modeling of regulatory elements is not stable yet.
           This means that this method might change or even disappear in a
           future release. Be aware of this if you use it.

 Returns : An array of Bio::SeqFeatureI implementing objects.
 Args    : 


=cut

sub promoters {
    my ($self) = @_;
    my @transcripts = $self->transcripts();
    my @feas = ();

    foreach my $tr (@transcripts) {
	push(@feas, $tr->promoters());
    }
    return @feas;
}


=head2 exons

 Title   : exons()
 Usage   : @exons = $gene->exons();
           @inital_exons = $gene->exons('Initial');
 Function: Get all exon features or all exons of a specified type of this gene
           structure.

           Exon type is treated as a case-insensitive regular expression and 
           optional. For consistency, use only the following types: 
           initial, internal, terminal, utr, utr5prime, and utr3prime. 
           A special and virtual type is 'coding', which refers to all types
           except utr.

           This method basically merges the exons returned by transcripts.

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects.
 Args    : An optional string specifying the type of exon.


=cut

sub exons {
    my ($self, @args) = @_;
    my @transcripts = $self->transcripts();
    my @feas = ();

    foreach my $tr (@transcripts) {
	push(@feas, $tr->exons(@args));
    }
    return @feas;
}

=head2 introns

 Title   : introns()
 Usage   : @introns = $gene->introns();
 Function: Get all introns of this gene structure.

           Note that this class currently generates these features on-the-fly,
           that is, it simply treats all regions between exons as introns.
           It assumes that the exons in the transcripts do not overlap.

           This method basically merges the introns returned by transcripts.

 Returns : An array of Bio::SeqFeatureI implementing objects.
 Args    : 


=cut

sub introns {
    my ($self) = @_;
    my @transcripts = $self->transcripts();
    my @feas = ();

    foreach my $tr (@transcripts) {
	push(@feas, $tr->introns());
    }
    return @feas;
}

=head2 poly_A_sites

 Title   : poly_A_sites()
 Usage   : @polyAsites = $gene->poly_A_sites();
 Function: Get the poly-adenylation sites of this gene structure.

           This method basically merges the poly-adenylation sites returned by
           transcripts.

 Returns : An array of Bio::SeqFeatureI implementing objects.
 Args    : 


=cut

sub poly_A_sites {
    my ($self) = @_;
    my @transcripts = $self->transcripts();
    my @feas = ();

    foreach my $tr (@transcripts) {
	push(@feas, $tr->poly_A_site());
    }
    return @feas;
}

=head2 utrs

 Title   : utrs()
 Usage   : @utr_sites = $gene->utrs('3prime');
           @utr_sites = $gene->utrs('5prime');
           @utr_sites = $gene->utrs();
 Function: Get the features representing untranslated regions (UTR) of this
           gene structure.

           You may provide an argument specifying the type of UTR. Currently
           the following types are recognized: 5prime 3prime for UTR on the
           5' and 3' end of the CDS, respectively.

           This method basically merges the UTRs returned by transcripts.

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects
           representing the UTR regions or sites.
 Args    : Optionally, either 3prime, or 5prime for the the type of UTR
           feature.


=cut

sub utrs {
    my ($self,@args) = @_;
    my @transcripts = $self->transcripts();
    my @feas = ();

    foreach my $tr (@transcripts) {
	push(@feas, $tr->utrs(@args));
    }
    return @feas;
}

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feats = $gene->sub_SeqFeature();
 Function: Returns an array of all subfeatures.

           This method is defined in Bio::SeqFeatureI. We override this here
           to include the transcripts.

 Returns : An array Bio::SeqFeatureI implementing objects.
 Args    : none


=cut

sub sub_SeqFeature {
    my ($self) = @_;   
    my @feas = ();

    # get what the parent already has
    @feas = $self->SUPER::sub_SeqFeature();
    push(@feas, $self->transcripts());
    return @feas;
}

=head2 flush_sub_SeqFeature

 Title   : flush_sub_SeqFeature
 Usage   : $gene->flush_sub_SeqFeature();
           $gene->flush_sub_SeqFeature(1);
 Function: Removes all subfeatures.

           This method is overridden from Bio::SeqFeature::Generic to flush
           all additional subfeatures, i.e., transcripts, which is
           almost certainly not what you want. To remove only features added
           through $gene->add_sub_SeqFeature($feature) pass any
           argument evaluating to TRUE.

 Example :
 Returns : none
 Args    : Optionally, an argument evaluating to TRUE will suppress flushing
           of all gene structure-specific subfeatures (transcripts).


=cut

sub flush_sub_SeqFeature {
   my ($self,$fea_only) = @_;

   $self->SUPER::flush_sub_SeqFeature();
   if(! $fea_only) {
       $self->flush_transcripts();
   }
}

sub gene_cleanup {
    my $self = shift;
    $self->flush_transcripts;
}

1;




