
#
# BioPerl module for Bio::SeqFeature::GeneStructure
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::GeneStructure

=head1 SYNOPSIS


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

Email hlapp@gmx.net or hilmar.lapp@pharma.novartis.com

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::GeneStructure;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::SeqFeature::Generic;


@ISA = qw(Bio::SeqFeature::Generic);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;
    
    my $make = $self->SUPER::_initialize(@args);

    my ($primary) =
	$self->_rearrange([qw(PRIMARY
			      )],@args);

    $primary = 'gene_structure' unless $primary;
    $self->primary_tag($primary);
    $self->strand(0) if(! defined($self->strand()));
    # set stuff in self from @args
    return $make; # success - we hope!
}

#
# Everything else is just inherited from SeqFeature::Generic.
#

=head2 promotors

 Title   : promotors()
 Usage   : @prom_sites = $gene->promotors();
 Function: Get the promotor features/sites of this gene structure.
 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           promotor regions or sites.
 Args    : 


=cut

sub promotors {
    my ($self) = @_;

    return undef unless exists($self->{'_promotors'});
    return @{$self->{'_promotors'}};
}

=head2 add_promotor

 Title   : add_promotor()
 Usage   : $gene->add_promotor($feature);
 Function: Add a promotor feature/site to this gene structure.
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.


=cut

sub add_promotor {
    my ($self, $fea) = @_;

    if(! $fea->isa('Bio::SeqFeatureI') ) {
	$self->throw("$fea does not implement Bio::SeqFeatureI");
    }
    if(! exists($self->{'_promotors'})) {
	$self->{'_promotors'} = [];
    }
    $self->_expand_region($fea);
    push(@{$self->{'_promotors'}}, $fea);
}

=head2 flush_promotors

 Title   : flush_promotors()
 Usage   : $gene->flush_promotors();
 Function: Remove all promotor features/sites from this gene structure.
 Returns : 
 Args    :


=cut

sub flush_promotors {
    my ($self) = @_;

    if(exists($self->{'_promotors'})) {
	delete($self->{'_promotors'});
    }
}

=head2 exons

 Title   : exons()
 Usage   : @exons = $gene->exons();
           @inital = $gene->exons('Initial');
 Function: Get all exon features or all exons of specified type of this gene
           structure.

           The optionally specified type simply screens the primary_tag of
           the exon features as the string it must begin with ('Initial'
           would match 'InitialExon' but not 'Exon'). Refer to the
           documentation of the class that produced this gene structure object
           for information about the possible types.

 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           exon regions.
 Args    : An optional string specifying the primary_tag of the feature.


=cut

sub exons {
    my ($self, $type) = @_;
    my @exons = ();

    return undef unless exists($self->{'_exons'});
    return @{$self->{'_exons'}} unless defined($type);
    foreach my $exon (@{$self->{'_exons'}}) {
	push(@exons, $exon) if($exon->primary_tag() =~ /^$type/);
    }
    return @exons;
}

=head2 add_exon

 Title   : add_exon()
 Usage   : $gene->add_exon($feature);
 Function: Add a exon feature to this gene structure.
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.


=cut

sub add_exon {
    my ($self, $fea, $type) = @_;
    my $key = '_exons';

    if(! $fea->isa('Bio::SeqFeatureI') ) {
	$self->throw("$fea does not implement Bio::SeqFeatureI");
    }
    if(! exists($self->{$key})) {
	$self->{$key} = [];
    }
    $self->_expand_region($fea);
    push(@{$self->{$key}}, $fea);
}

=head2 flush_exons

 Title   : flush_exons()
 Usage   : $gene->flush_exons();
 Function: Remove all exon features from this gene structure.
 Returns : 
 Args    :


=cut

sub flush_exons {
    my ($self) = @_;

    if(exists($self->{'_exons'})) {
	delete($self->{'_exons'});
    }
}

=head2 introns

 Title   : introns()
 Usage   : @introns = $gene->introns();
 Function: Get all intron features this gene structure.

           Note that this class currently generates these features on-the-fly,
           that is, it simply treats all regions between exons as introns.
           It assumes that the exons were added consecutively and that they
           do not overlap. Hence, it is merely a convenience method, but this
           may change in the future or in overriding classes.

 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           intron regions.
 Args    : 


=cut

sub introns {
    my ($self) = @_;
    my @introns = ();
    my @exons = $self->exons();

    return undef unless @exons;
    for(my $i = 0; $i < $#exons; $i++) {
	my ($start, $end);
	my $intron;

	if($exons[$i]->end() < $exons[$i+1]->start()) {
	    # normal order
	    $start = $exons[$i]->end() + 1;
	    $end = $exons[$i+1]->start() - 1;
	} else {
	    # reversed order
	    $start = $exons[$i+1]->end() + 1;
	    $end = $exons[$i]->start() - 1;
	}
	$intron = Bio::SeqFeature::Generic->new(
                                        '-start'   => $start,
                                        '-end'     => $end,
                                        '-strand'  => $exons[$i]->strand(),
                                        '-primary' => 'Intron',
					'-source'  => 'intervening');
	my $seq = $self->entire_seq();
	$intron->attach_seq($seq) if $seq;
	push(@introns, $intron);
    }
    return @introns;
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

    return undef unless exists($self->{'_poly_A_sites'});
    return @{$self->{'_poly_A_sites'}};
}

=head2 add_poly_A_site

 Title   : add_poly_A_site()
 Usage   : $gene->add_poly_A_site($feature);
 Function: Add a poly-adenylation feature/site to this gene structure.
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.


=cut

sub add_poly_A_site {
    my ($self, $fea) = @_;

    if(! $fea->isa('Bio::SeqFeatureI') ) {
	$self->throw("$fea does not implement Bio::SeqFeatureI");
    }
    if(! exists($self->{'_poly_A_sites'})) {
	$self->{'_poly_A_sites'} = [];
    }
    $self->_expand_region($fea);
    push(@{$self->{'_poly_A_sites'}}, $fea);
}

=head2 flush_poly_A_sites

 Title   : flush_poly_A_sites()
 Usage   : $gene->flush_poly_A_sites();
 Function: Remove all poly-adenylation features/sites from this gene structure.
 Returns : 
 Args    :


=cut

sub flush_poly_A_sites {
    my ($self) = @_;

    if(exists($self->{'_poly_A_sites'})) {
	delete($self->{'_poly_A_sites'});
    }
}

=head2 utrs

 Title   : utrs()
 Usage   : @utr_sites = $gene->utrs();
 Function: Get the UTR features/sites of this gene structure.
 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           UTR regions or sites.
 Args    : 


=cut

sub utrs {
    my ($self) = @_;

    return undef unless exists($self->{'_utrs'});
    return @{$self->{'_utrs'}};
}

=head2 add_utr

 Title   : add_utr()
 Usage   : $gene->add_utr($feature);
 Function: Add a UTR feature/site to this gene structure.
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.


=cut

sub add_utr {
    my ($self, $fea) = @_;

    if(! $fea->isa('Bio::SeqFeatureI') ) {
	$self->throw("$fea does not implement Bio::SeqFeatureI");
    }
    if(! exists($self->{'_utrs'})) {
	$self->{'_utrs'} = [];
    }
    $self->_expand_region($fea);
    push(@{$self->{'_utrs'}}, $fea);
}

=head2 flush_utrs

 Title   : flush_utrs()
 Usage   : $gene->flush_utrs();
 Function: Remove all UTR features/sites from this gene structure.
 Returns : 
 Args    :


=cut

sub flush_utrs {
    my ($self) = @_;

    if(exists($self->{'_utrs'})) {
	delete($self->{'_utrs'});
    }
}

=head2 seqname

 Title   : seqname
 Usage   : $gene->seqname($seq)
 Function: Sets the name of the sequence this feature applies to.

           We override this here in order to ensure that the name is set
           for all subfeatures like exons etc.

 Example :
 Returns : 
 Args    :


=cut

sub seqname {
    my ($self, $val) = @_;
    my $fea;

    # first let the parent do its work
    my $ret = $self->SUPER::seqname($val);

    if(defined($val)) {
	# set for subfeaturess if they want
	foreach $fea ($self->exons()) {
	    if($fea->can("seqname") ) {
		$fea->seqname($val);
	    }
	}
	foreach $fea ($self->promotors()) {
	    if($fea->can("seqname") ) {
		$fea->seqname($val);
	    }
	}
	foreach $fea ($self->poly_A_sites()) {
	    if($fea->can("seqname") ) {
		$fea->seqname($val);
	    }
	}
	foreach $fea ($self->utrs()) {
	    if($fea->can("seqname") ) {
		$fea->seqname($val);
	    }
	}
    }
    return $ret;
}

=head2 attach_seq

 Title   : attach_seq
 Usage   : $gene->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence.

           We override this here in order to ensure that the seq is attached
           to all subfeatures like exons etc.

 Example :
 Returns : TRUE on success
 Args    :


=cut

sub attach_seq {
    my ($self, $seq) = @_;
    my $fea;

    # first let the parent do its work
    my $ret = $self->SUPER::attach_seq($seq);

    # attach to subfeatures if they want it
    foreach $fea ($self->exons()) {
	if($fea->can("attach_seq") ) {
	    $fea->attach_seq($seq);
	}
    }
    foreach $fea ($self->promotors()) {
	if($fea->can("attach_seq") ) {
	    $fea->attach_seq($seq);
	}
    }
    foreach $fea ($self->poly_A_sites()) {
	if($fea->can("attach_seq") ) {
	    $fea->attach_seq($seq);
	}
    }
    foreach $fea ($self->utrs()) {
	if($fea->can("attach_seq") ) {
	    $fea->attach_seq($seq);
	}
    }
    return $ret;
}

=head2 cds

 Title   : cds
 Usage   : $cds_dna = $gene->cds();
 Function: Returns the CDS (coding sequence) as defined by this gene structure
           and the attached sequence.

           If no sequence is attached this method will throw an exception.

           Note that the implementation provided here returns a concatenation
           of all exons, thereby assuming that these were added consecutively
           and that exons do not overlap. Hence, it is merely a convenience
           method, but this may be changed by overriding classes/methods.

           Note also that you cannot set the CDS via this method. It will
           throw an exception on an attempt to do so. Set a single CDS
           feature as an exon, or derive your own class if you want to
           store a predicted CDS.

 Example :
 Returns : A string holding the DNA sequence defined as coding by the regions
           of the exons.
 Args    :

=cut

sub cds {
    my ($self, @args) = @_;
    my @exons = $self->exons();
    my $cds = "";

    if(! $self->entire_seq()) {
	$self->throw("CDS requested, but no sequence object attached");
    }
    if(@args) {
	$self->throw("Unexpected argument. You cannot set the CDS via this method");
    }
    if($#exons > 0) {
	# determine order
	if((($exons[0]->strand() == -1) &&
	    ($exons[0]->end() < $exons[1]->start())) ||
	   (($exons[0]->strand() == 1) &&
	    ($exons[0]->start() > $exons[1]->end()))) {
	    # normal order, but minus strand, or reverse order, but plus
	    # strand: we need to reverse the order in either case
	    @exons = CORE::reverse(@exons);
	}
	foreach my $exon (@exons) {
	    $cds .= $exon->seq();
	}
    }
    return $cds;
}

sub _tag_value {
    my ($self, $tag, $value) = @_;

    if(defined($value)) {
	$self->remove_tag($tag) if($self->has_tag($tag));
	$self->add_tag_value($tag, $value);
    }
    return ($self->each_tag_value($tag))[0];
}

1;
