
#
# BioPerl module for Bio::AnnSeq
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnSeq - Annotated Sequence

=head1 SYNOPSIS

    $stream = Bio::AnnSeqIO->new(-file 'my.embl',-format => 'EMBL')

    foreach $annseq ( $stream->next_annseq() ) {
	foreach $feat ( $annseq->all_SeqFeatures() ) {
	    print "Feature ",$feat->primary_tag," at ", $feat->start, " ",$feat->end, "\n";
	}
    }

=head1 DESCRIPTION

An AnnSeq is a sequence with sequence features placed on them

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnnSeq;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::Annotation;
use Bio::Seq;

@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  my($ann);
  my $make = $self->SUPER::_initialize;
  $self->{'_as_feat'} = [];
  $ann = new Bio::Annotation;
  $self->annotation($ann);

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: 
 Example : 
 Returns : value of seq
 Args    : newvalue (optional)


=cut

sub seq{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'seq'} = $value;
      # descend down over all seqfeature objects, seeing whether they
      # want an attached seq.

      foreach my $sf ( $obj->top_SeqFeatures() ) {
	  if( $sf->can("attach_seq") ) {
	      $sf->attach_seq($value);
	  }
      }

    }
    return $obj->{'seq'};

}

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($seq_obj)
 Function: 
 Example : 
 Returns : value of annotation
 Args    : newvalue (optional)


=cut

sub annotation{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

}

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $annseq->add_SeqFeature($feat);
 Function: Adds t
 Example :
 Returns : 
 Args    :


=cut

sub add_SeqFeature{
   my ($self,$feat) = @_;
   my ($fseq,$aseq);

   if( !$feat->isa("Bio::SeqFeatureI") ) {
       $self->warn("$feat is not a SeqFeatureI and that's what we expect...");
   }

   if( $feat->can("seq") ) {
       $fseq = $feat->seq;
       $aseq = $self->seq;

       if( defined $aseq ) {
	   if( defined $fseq ) {
	       if( $aseq ne $fseq ) {
		   $self->warn("$feat has an attached sequence which is not in this annseq. I worry about this");
	       }
	   } else {
	       if( $feat->can("attach_seq") ) {
		   # attach it 
		   $feat->attach_seq($aseq);
	       }
	   }
       } # end of if aseq
   } # end of if the feat can seq

   push(@{$self->{'_as_feat'}},$feat);
       
}

=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub top_SeqFeatures{
   my ($self) = @_;

   return @{$self->{'_as_feat'}};
}

=head2 all_SeqFeatures

 Title   : all_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub all_SeqFeatures{
   my ($self) = @_;
   my (@array);
   foreach my $feat ( $self->top_SeqFeatures() ){
       push(@array,$feat);
       &_retrieve_subSeqFeature(\@array,$feat);
   }

   return @array;
}


sub _retrieve_subSeqFeature {
    my ($arrayref,$feat) = @_;

    foreach my $sub ( $feat->sub_SeqFeature() ) {
	push(@$arrayref,$sub);
	&_retrieve_subSeqFeature($arrayref,$sub);
    }

}

=head2 fetch_SeqFeatures

 Title   : fetch_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_SeqFeatures{
   my ($self,@args) = @_;

   $self->throw("Not implemented yet");
}


=head2 species

 Title   : species
 Usage   : 
 Function: Gets or sets the species
 Example : $species = $self->species();
 Returns : Bio::Species object
 Args    : Bio::Species object or none;


=cut

sub species {
    my ($self, $species) = @_;

    if ($species) {
        $self->{'species'} = $species;
    } else {
        return $self->{'species'}
    }
}









