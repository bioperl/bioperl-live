
#
# BioPerl module for Bio::Pfam::Annotation::Comment
#
# Cared for by Ewan Birney <pfam@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Comment - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::Comment;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'flat'} = [];
# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 text

 Title   : text
 Usage   : $self->text($newval)
 Function: 
 Example : 
 Returns : value of text
 Args    : newvalue (optional)


=cut

sub text{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'text'} = $value;
    }
    return $self->{'text'};

}



1;
