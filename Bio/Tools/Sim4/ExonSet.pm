
#
# BioPerl module for Bio::Tools::Sim4::ExonSet
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Sim4::ExonSet - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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


package Bio::Tools::Sim4::ExonSet;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::SeqFeature::Generic;

@ISA = qw(Bio::SeqFeature::Generic);

# new() is inherited from Bio::Root::Object

# _initial2ize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'exon_array'} = [];
  $self->primary_tag('ExonSet');
  $self->source_tag('Sim4');
  $self->strand(0);
# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 sub_SeqFeature
    
 Title   : sub_SeqFeature
 Usage   : foreach $feat ( $exonset->sub_SeqFeature )
 Function:
 Example :
 Returns : 
 Args    :

  Note: we overload this to add in the Exons stored
  separately from each_Exon, but also allow people to
  add their own sub_SeqFeatures

=cut

sub sub_SeqFeature{
   my ($self) = @_;

   my @feat = $self->SUPER::sub_SeqFeature();

   push(@feat,$self->each_Exon());

   return @feat;
}

=head2 each_Exon

 Title   : each_Exon
 Usage   : foreach $exon ( $exonset->each_Exon ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Exon{
   my ($self) = @_;


   return @{$self->{'exon_array'}};
}

=head2 add_Exon

 Title   : add_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Exon{
   my ($self,$exon) = @_;

   $exon->isa('Bio::Tools::Sim4::Exon') || $self->throw("$exon is not a Bio::Tools::Sim4::Exon. Yuk!");

   # expand the start/end point as necessary.

   if( !defined $self->start && !defined $self->end ) {
       $self->start($exon->start());
       $self->end($exon->end());
       $self->strand($exon->strand);
   } else {
       my ($start,$end,$strand) = $self->union($exon);
       $self->start($start);
       $self->end($end);
       $self->strand($strand);
   }

   push(@{$self->{'exon_array'}},$exon);
}


1;
