# $Id$
#
# BioPerl module for Bio::SearchIO::EventHandlerI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::EventHandlerI - An abstract Event Handler

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


package Bio::SearchIO::EventHandlerI;
use vars qw(@ISA);
use strict;
use Carp;

use Bio::Event::EventHandlerI;
use Bio::Event::EventGeneratorI;

@ISA = qw (Bio::Event::EventHandlerI);

=head2 start_report

 Title   : start_report
 Usage   : $handler->start_report($data)
 Function: Begins a report event cycle
 Returns : none 
 Args    : Type of Report

=cut

sub start_report {
    my ($self) = @_;
    $self->_abstractDeath('start_report');
}

=head2 end_report

 Title   : end_report
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_report{
   my ($self,@args) = @_;
   $self->_abstractDeath('end_report');
}

=head2 start_hsp

 Title   : start_hsp
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_hsp{
   my ($self,@args) = @_;
   $self->_abstractDeath('start_hsp');
}

=head2 end_hsp

 Title   : end_hsp
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_hsp{
   my ($self,@args) = @_;
   $self->_abstractDeath('end_hsp');
}

=head2 start_subject

 Title   : start_subject
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_subject{
   my ($self,@args) = @_;
   $self->_abstractDeath('start_subject');
}

=head2 end_subject

 Title   : end_subject
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_subject{
   my ($self,@args) = @_;
   $self->_abstractDeath('end_subject');
}

=head2 Bio::Event::EventHandlerI methods

=head2 will_handle

 Title   : will_handle
 Usage   : if( $handler->will_handle($event_type) ) { ... }
 Function: Tests if this event builder knows how to process a specific event
 Returns : boolean
 Args    : event type name


=cut

=head2 SAX methods

=head2 start_document

 Title   : start_document
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=head2 end_document

 Title   : end_document
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=head2 start_element

 Title   : start_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=head2 end_element

 Title   : end_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=head2 in_element

 Title   : in_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=head2 within_element

 Title   : within_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=head2 characters

 Title   : characters
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

1;
