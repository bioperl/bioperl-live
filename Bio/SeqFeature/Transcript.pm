
#
# BioPerl module for Bio::SeqFeature::Transcript
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Transcript - A particular transcription/processing of a gene

=head1 SYNOPSIS

   # PLEASE READ Bio::SeqFeature::Gene documentation
   # before trying to decode this. The gene object is
   # the key object to use to manipulate genes etc

   $trans = Bio::Transcript->new();
   $trans->add_Exon(10,20,-1); # start,end,strand
   # add more exons. Exons can be added in any order.

   # creates a translation from this transcript at this particular
   # start/stop.

   $tls = $trans->create_Translation(140,900); # start/end in genomic coordinates
   $gene->add_Translation($tls);
    
   # implicitly indicates that this is a protein coding gene
   $trans->Translation($tls);

   

   foreach $trans ( $gene->each_Transcript() ) {
       foreach $exon ( $trans->each_Exon() ) {
	   # do something with exons
       }
       $cdna = $trans->seq();
       # not all transcripts necessarily have translations
       if( $trans->is_protein_coding ) {
	   $pep = $trans->Translation->seq();
       }
   }


=head1 DESCRIPTION


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


package Bio::SeqFeature::Transcript;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Translation;
use Bio::SeqFeature::Exon;


@ISA = qw(Bio::SeqFeature::Generic);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'_trans_exon_ary'} = [];
# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 Translation

 Title   : Translation
 Usage   : $obj->Translation($newval)
           $tls = $obj->Translation()
 Function: 
 Example : 
 Returns : value of Translation
 Args    : newvalue (optional)


=cut

sub Translation{
   my ($obj,$value) = @_;

   if( defined $value) {
       if( ! ref $value || ! $value->isa('Bio::SeqFeature::Translation') ) {
	   $obj->throw("set_Translation requires one argument being the translation, not $value");
       }
       $obj->{'_trans_translation'} = $value;
   }
   return $obj->{'_trans_translation'};
}

=head2 is_protein_coding

 Title   : is_protein_coding
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub is_protein_coding{
   my ($self) = @_;

   return exists $self->{'_trans_translation'};
}

=head2 each_Exon

 Title   : each_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Exon{
   my ($self) = @_;

   return @{$self->{'_trans_exon_ary'}};
}

=head2 add_Exon

 Title   : add_Exon
 Usage   : $trans->add_Exon($exon);
           $trans->add_Exon(@exons);
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Exon{
   my ($self) = shift @_;

   foreach my $exon ( @_ ) {
       if( ! ref $exon ) {
	   $self->throw("add_Exon(\$exon) assummes that exon is an object");
       }
       if( ! $exon->isa('Bio::SeqFeatureI') ) {
	   $self->warn("Exons have to be at least SeqFeatureI's - which should be easy to add to $exon");
       }
       push(@{$self->{'_trans_exon_ary'}},$exon);
       
   }

   $self->_sort_Exons();
}


=head2 seq

 Title   : seq
 Usage   : $cdna = $trans->seq();
 Function: Returns a Bio::Seq object of the cDNA
 Example :
 Returns : 
 Args    :


=cut

sub seq{
   my ($self) = @_;
   my $str;
   foreach my $exon ( $self->each_Exon() ) {
       my $tseq = $exon->seq();
       $str .= $tseq->str();
   }

   # FIXME - should not be cdna here!
   my $out = Bio::Seq->new(-seq => $str , -id => "cdna" );
}
=head2 _attach_seq

 Title   : _attach_seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _attach_seq{
   my ($self,$seq) = @_;

   foreach my $exon ( $self->each_Exon() ) {
       $exon->attach_seq($seq);
   }
}

=head2 _sort_Exons

 Title   : _sort_Exons
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _sort_Exons{
   my ($self) = @_;

   @{$self->{'_trans_exon_ary'}} = sort { $a->start() <=> $b->start() } @{$self->{'_trans_exon_ary'}};
}

=head2 create_Translation

 Title   : create_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub create_Translation{
   my ($self,$start,$end) = @_;
   my $tls;
   if( !defined $end ) {
       $self->throw("create_Translation requires start/end points");
   }

   $tls = Bio::SeqFeature::Translation->new();
   # find first exon
   my $first;
   foreach my $exon ( $self->each_Exon() ) {
       if( $exon->start <= $start && $exon->end() >= $start ) {
	   $first = $exon;
	   last;
       }
   }

   if( ! defined $first ) {
       $self->throw("$start is not in the exons of the transcript (remember the positions have to be in absolute coordinates)");
   }

   # find last exon
   my $last;
   foreach my $exon ( $self->each_Exon() ) {
       if( $exon->start <= $end && $exon->end() >= $end ) {
	   $last = $exon;
	   last;
       }
   }

   if( ! defined $last) {
       $self->throw("$end is not in the exons of the transcript (remember the positions have to be in absolute coordinates)");
   }

   $tls->first_exon(Bio::SeqFeature::Exon->new(-start => $start, -end => $first->end, -strand => $first->strand));
   $tls->last_exon(Bio::SeqFeature::Exon->new(-start => $last->start, -end => $end, -strand => $last->strand));
   
   foreach my $exon ( $self->each_Exon() ) {
       if( $exon->start > $start && $exon->end() < $end ) {
	   $tls->_add_internal_exon($exon);
       }
   }
		   
   return $tls;
}




1;








