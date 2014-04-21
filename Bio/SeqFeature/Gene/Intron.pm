#
# BioPerl module for Bio::SeqFeature::Gene::Intron
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

Bio::SeqFeature::Gene::Intron - An intron feature

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene::Intron;
use strict;

use Bio::SeqFeature::Gene::Exon;

use base qw(Bio::SeqFeature::Gene::NC_Feature);

sub new {
    my($class,@args) = @_;

    # introns are non-coding by default
    if(! grep { lc($_) eq '-is_coding'; } @args) {
	push(@args, '-is_coding', 0);
    }
    my $self = $class->SUPER::new(@args);

    my ($primary, $prim) = 
	$self->_rearrange([qw(PRIMARY PRIMARY_TAG)],@args);
    $self->primary_tag('intron') unless $primary || $prim;

    return $self;
}

=head2 upstream_Exon 

  Title   : upstream_Exon 
  Usage   : $intron->upstream_Exon()
  Function: exon upstream of the intron 
  Returns : Bio::EnsEMBL::Exon
  Args    : 

=cut

sub upstream_Exon {
    my( $self, $exon ) = @_;
    
    if ($exon) {
        $self->{'_intron_location'} = undef;
        $self->throw("'$exon' is not a Bio::SeqFeature::Gene::ExonI") 
	    unless $exon->isa('Bio::SeqFeature::Gene::ExonI');
        $self->{'_upstream_exon'} = $exon;
    }
    return $self->{'_upstream_exon'};
}


=head2 downstream_Exon 

  Title   : downstream_Exon 
  Usage   : $intron->downstream_Exon()
  Function: exon downstream of the intron 
  Returns : Bio::EnsEMBL::Exon
  Args    : 

=cut

sub downstream_Exon {
    my( $self, $exon ) = @_;
    
    if ($exon) {
        $self->{'_intron_location'} = undef;
        $self->throw("'$exon' is not a Bio::SeqFeature::Gene::ExonI")
            unless $exon->isa('Bio::SeqFeature::Gene::ExonI');
        $self->{'_downstream_exon'} = $exon;
    }
    return $self->{'_downstream_exon'};
}

=head2 phase 

  Title   : phase 
  Usage   : $intron->phase()
  Function: returns the phase of the intron(where it interrupts the codon)  
  Returns : int(0,1,2)
  Args    : 

=cut

sub phase {
  my ($self) = @_;
  return $self->downstream_Exon->phase;
}


=head2 acceptor_splice_site 

  Title   : acceptor_splice_site 
  Usage   : $intron->acceptor_splice_site(21,3)
  Function: returns the sequence corresponding to the 
            consensus acceptor splice site. If start and
            end are provided, it will number of base pairs
            left and right of the canonical AG. Here 21 means
            21 bp into intron and 3 means 3 bp into the exon.
            --Intron--21----|AG|-3-----Exon
            Defaults to 21,3

  Returns : Bio::Seq
  Args    : start and end

=cut

sub acceptor_splice_site {
  my ($self,$ss_start,$ss_end) = @_;
  $ss_start = 21 unless defined $ss_start;
  $ss_end   = 3 unless defined $ss_end;
  if($self->strand < 0){
    my $tmp= $ss_start;
    $ss_start = $ss_end;
    $ss_end = $tmp;
  }
  my $intron_end= $self->location->end;
  my $down_exon = $self->downstream_Exon;
  my $acceptor;  
  if($self->strand < 0){
      $ss_start= $ss_start >  $down_exon->length ? $down_exon->length: $ss_start;
      $ss_end= $ss_end > $self->length-2 ? $self->length-2 : $ss_end;
      $acceptor = Bio::SeqFeature::Generic->new(-start=>$self->start - ($ss_start) ,  
                                               -end=>$self->start + ($ss_end+1),
                                               -strand=>$self->strand,    
                                               -primary_tag=>"donor splice site");
  } 
  else {
    $ss_start = $ss_start > $self->length-2 ? $self->length-2 : $ss_start;
    $ss_end = $ss_end > $down_exon->length ? $down_exon->length : $ss_end;
 

    $acceptor = Bio::SeqFeature::Generic->new(-start=>$self->end - ($ss_start + 1),  
                                               -end=>$self->end + $ss_end,
                                               -strand=>$self->strand,    
                                               -primary_tag=>"donor splice site");
  }
  $acceptor->attach_seq($self->entire_seq);

  return $acceptor;

}

=head2 donor_splice_site 

  Title   : donor_splice_site 
  Usage   : $intron->donor_splice_site(3,6)
  Function: returns the sequence corresponding to the 
            consensus donor splice site. If start and
            end are provided, it will number of base pairs
            left and right of the canonical GT. Here 3 means
            3 bp into exon and 6 means 6 bp into the intron.
            --Exon-3--|GT|-6----Intron-
            Defaults to 3,6

  Returns : Bio::Seq
  Args    : start and end

=cut

sub donor_splice_site {
  my ($self,$ss_start,$ss_end) = @_;
  $ss_start = 3 unless defined $ss_start;
  $ss_end   = 10 unless defined $ss_end;
  if($self->strand < 0){
    my $tmp= $ss_start;
    $ss_start = $ss_end;
    $ss_end = $tmp;
  }
  my $up_exon = $self->upstream_Exon;
  my $donor;
  if($self->strand < 0){
    $ss_end = $ss_end > $up_exon->length ? $up_exon->length : $ss_end;
    $ss_start   = $ss_start> $self->length -2 ? $self->length -2 : $ss_start;
    $donor = Bio::SeqFeature::Generic->new(-start=>$self->end -  ($ss_start+1),
                                            -end  => $self->end + ($ss_end),
                                            -strand=>$self->strand,
                                            -primary_tag=>"acceptor splice site");
  } 
  else {
    $ss_start = $ss_start > $up_exon->length ? $up_exon->length : $ss_start;
    $ss_end   = $ss_end > $self->length -2 ? $self->length -2 : $ss_end;
    $donor = Bio::SeqFeature::Generic->new(-start=>$self->start - $ss_start,
                                           -end  => $self->start +($ss_end+1),
                                            -strand=>$self->strand,
                                            -primary_tag=>"acceptor splice site");
  }
  $donor->attach_seq($self->entire_seq);
  return $donor;
} 

sub location {
    my( $self ) = @_;
    
    unless ($self->{'_intron_location'}) {
        my $loc = Bio::Location::Simple->new;
    
        my $up_exon = $self->upstream_Exon;
        my $down_exon = $self->downstream_Exon;
        
        # Get the PrimarySeqs attached to both and check it is the same sequence
        my $up_seq   = $up_exon  ->entire_seq;
        my $down_seq = $down_exon->entire_seq;
        unless (ref($up_seq) eq ref($down_seq) ) {
            $self->throw("upstream and downstream exons are attached to different sequences\n'$up_seq' and '$down_seq'");
        }
        
        # Check that the exons are on the same strand.  (Do I need to bother?)
        my $up_strand   = $up_exon  ->strand;
        my $down_strand = $down_exon->strand;
        unless ($up_strand == $down_strand) {
            $self->throw("upstream and downstream exons are on different strands "
                . "('$up_strand' and '$down_strand')");
        }
        $loc->strand($up_strand);
        
        #   $exon_end is the  end  of the exon which is 5' of the intron on the genomic sequence.
        # $exon_start is the start of the exon which is 3' of the intron on the genomic sequence.
        my( $exon_end, $exon_start );
        if ($up_strand == 1) {
            $exon_end   = $up_exon  ->end;
            $exon_start = $down_exon->start;
        } else {
            $exon_end   = $down_exon->end;
            $exon_start = $up_exon  ->start;
        }
        unless ($exon_end < $exon_start) {
            $self->throw("Intron gap begins after '$exon_end' and ends before '$exon_start'");
        }
        $loc->start($exon_end   + 1);
        $loc->end  ($exon_start - 1);
        
        # Attach the sequence and location objects to the intron
        $self->{'_intron_location'} = $loc;
        
    }
    return $self->{'_intron_location'};
}
1;
