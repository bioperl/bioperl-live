
#
# BioPerl module for Bio::SeqFeature::Homol
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Homol - Sequence Feature with Homolgous tags

=head1 SYNOPSIS

    $sf = new Bio::SeqFeature::Homol ( -start => 10, 
				       -end   => 50,
				       -strand => 1,
				       );
    $sf->attach_seq($seq1);

    $hsf = new Bio::SeqFeature ( -start => 405,
				 -end => 410,
				 -strand => 1,
				 );

    $hsf->attach_seq($seq2);



    $sf->homol_SeqFeature($hsf);

=head1 DESCRIPTION

This is a derived object off SeqFeature::Generic, which provides a homol slot, allowing
this feature to be linked to a SeqFeature presumably on a different sequence with the
concept that the two seqfeatures are somehow homologous. 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Homol;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::SeqFeature::Generic;


@ISA = qw(Bio::SeqFeature::Generic);

# new() is inherited from Bio::SeqFeature::Generic

sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize(@args);

  return $make;
}

=head2 homol_SeqFeature

 Title   : homol_SeqFeature
 Usage   : $obj->homol_SeqFeature($seqfeature)
 Function: Provides a set/get onto the homol link from this SeqFeature to
           another SeqFeature on a different sequence presummed to be homologous
           to this one. 
 Returns : homol_SeqFeature (if any)
 Args    : newvalue (optional)

  Notes:

    This does not work properely when $seqfeature is a Bio::SeqFeature::Homol itself,
  in which case we should put in the reciprocal link - but that gives us a nasty 
  circular reference.

  Also this should coordinate with Aarons stuff. I have to look at that...

=cut

sub homol_SeqFeature{
   my $self = shift;
   if( @_ ) {
       my $value = shift;
       $value->isa("Bio::SeqFeatureI") || $self->throw("$value is not a Bio::SeqFeatureI implementing object");
       $self->{'homol_SeqFeature'} = $value;
       if( $value->isa("Bio::SeqFeature::Homol") ) {
	   my $h = $value->homol_SeqFeature();
	   if( ! defined $h ) {
	       $value->homol_SeqFeature($self);
	   } else {
	       if( $h != $self ) {
		   $self->throw("Attempting to attach a homol seqfeature that already has a homol partner, but not myself!");
	       }
	   }

       } else {
	   $self->warn("You are now suggested to use Homol objects paired with Homol objects, rather than generic seqfeatures.");
       }
   }
   return $self->{'homol_SeqFeature'};

}

=head2 _remove_homol_SeqFeature

 Title   : _remove_homol_SeqFeature
 Usage   :
 Function: internal remove function, mainly to keep -w happy
           about removing things
 Example :
 Returns : 
 Args    :


=cut

sub _remove_homol_SeqFeature{
   my ($self) = @_;

   $self->{'homol_SeqFeature'} = undef;
}


=head2 DESTROY

 Title   : DESTROY
 Usage   : Object deconstructor. We need one to get rid of the potential 
           self reference in the homol pair

=cut

sub DESTROY {
    my $self = shift;
    my $h = $self->homol_SeqFeature();

    if( defined $h ) {
	if( $h->can('homol_SeqFeature') ) {
	    my $partner  = $h->homol_SeqFeature();
	    if( defined $partner && $partner == $self ) {
		$h->_remove_homol_SeqFeature();
	    }
	}
    }

    $self->SUPER::DESTROY();
}



