# $Id$
#
# BioPerl module for Bio::Event::EventGeneratorI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Event::EventGeneratorI - DESCRIPTION of Interface

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Event::EventGeneratorI;
use vars qw(@ISA);
use strict;
use Carp;

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  confess "Abstract method '$caller' defined in interface Bio::Event::EventGeneratorI not implemented by package $package. Not your fault - author of $package should be blamed!";
}


=head2 SAX methods

=head2 start_document

 Title   : start_document
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_document{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 end_document

 Title   : end_document
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_document{
   my ($self,@args) = @_;
   $self->_abstractDeath;

}

=head2 start_element

 Title   : start_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_element{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 end_element

 Title   : end_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_element{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}


=head2 in_element

 Title   : in_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub in_element{
   my ($self,@args) = @_;
   $self->_abstractDeath;

}

=head2 within_element

 Title   : within_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub within_element{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 characters

 Title   : characters
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub characters{
   my ($self,@args) = @_;
   $self->_abstractDeath;

}

1;
