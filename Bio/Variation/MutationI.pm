
#
# BioPerl module for Bio::Variation::MutationI
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Variation::MutationI - Sequence Change SeqFeature abstract class

=head1 SYNOPSIS

    if ($obj->isa('Bio::Variation::MutationI')) {
	do something;
    }

=head1 DESCRIPTION

This abstract superclass defines a simplified interface to be sequence
variations describing canonical mutations from the wild type default
allele into the mutated allele. The compliations of potentially
handling several alleles in a polymorphic site can be ignored.

Also, it disables the default
L<Bio::Variation::VariantI> methods add_Allele and each_allele.


=head1 CONTACT

Heikki Lehvaslaiho <heikki@ebi.ac.uk>

=cut

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Variation::MutationI;
use vars qw(@ISA);
use strict;

# Object preamble - inheritance

#use Bio::Variation::Allele;
#@ISA = qw( Bio::Variation::VariantI);


sub add_Allele {
  my ($self) = @_;
  $self->throw("Mutations are restrited to using methods allele_ori and allele_mut, only");
}

sub each_Allele{
  my ($self) = @_;
  $self->throw("Mutations are restrited to using methods allele_ori and allele_mut, only");
}


=head2 allele_ori

 Title   : allele_ori
 Usage   : $obj->allele_ori();
 Function: 

            Sets and returns original allele sequence.  If value is not
            set, returns false. 

            Amino acid sequences are stored in upper case characters,
            others in lower case.

 Example : 
 Returns : string
 Args    : string

=cut


sub allele_ori {
   my ($obj,$value) = @_;
   if( defined $value) {
       if ($value->isa('Bio::Variation::Allele')) {
	   $obj->throw("Value is not Bio::Variation::Allele but [$value]");
       } else {
	   $obj->{'allele_ori'} = $value;
       }
   }
   return $obj->{'allele_ori'};
}


=head2 allele_mut

 Title   : allele_mut
 Usage   : $obj->allele_mut();
 Function: 

            Sets and returns the mutated allele sequence.  If value is not
            set, returns false. 

 Example : 
 Returns : string
 Args    : string

=cut


sub allele_mut {
   my ($obj,$value) = @_;
   if( defined $value) {
       if ($value->isa('Bio::Variation::Allele')) {
	   $obj->throw("Value is not Bio::Variation::Allele but [$value]");
       } else {
	   $obj->{'allele_mut'} = $value;
       }
   }
   return $obj->{'allele_mut'};
}




1;








