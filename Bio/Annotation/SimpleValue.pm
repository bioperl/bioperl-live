
#
# BioPerl module for Bio::Annotation::SimpleValue
#
# Cared for by bioperl <bioperl-l@bio.perl.org>
#
# Copyright bioperl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::SimpleValue - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Scalar value annotation object 

=head1 AUTHOR - bioperl

Email bioperl-l@bio.perl.org

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::SimpleValue;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::AnnotationI;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI Bio::AnnotationI);
# new() can be inherited from Bio::Root::RootI

=head2 AnnotationI implementing functions

=cut

=head2 as_text

 Title   : as_text
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub as_text{
   my ($self) = @_;

   return "Value: ".$self->value;
}

=head2 hash_tree

 Title   : hash_tree
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub hash_tree{
   my ($self) = @_;
   
   my $h = {};
   $h->{'value'} = $self->value;
}

=head2 Specific accessors for SimpleValue

=cut

=head2 value

 Title   : value
 Usage   : $obj->value($newval)
 Function: 
 Example : 
 Returns : value of value
 Args    : newvalue (optional)


=cut

sub value{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'value'} = $value;
    }
    return $self->{'value'};

}

1;
