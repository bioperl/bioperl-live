#
# BioPerl module for Bio::Event::EventHandlerI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Event::EventHandlerI - An Event Handler Interface

=head1 SYNOPSIS

    # do not use this module directly
    # See Bio::SearchIO::SearchResultEventHandler for an example of
    # implementation.

=head1 DESCRIPTION

This interface describes the basic methods required for
EventHandlers.  These are essentially SAX methods. 

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

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Event::EventHandlerI;
use strict;
use Carp;

use base qw(Bio::Root::RootI);

=head2 will_handle

 Title   : will_handle
 Usage   : if( $handler->will_handle($event_type) ) { ... }
 Function: Tests if this event builder knows how to process a specific event
 Returns : boolean
 Args    : event type name


=cut

sub will_handle{
   my ($self,$type) = @_;
   $self->throw_not_implemented();
}

=head2 SAX methods

=cut

=head2 start_document

 Title   : start_document
 Usage   : $eventgenerator->start_document();
 Function: Handle a start document event
 Returns : none
 Args    : none


=cut

sub start_document{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 end_document

 Title   : end_document
 Usage   : $eventgenerator->end_document();
 Function: Handle an end document event
 Returns : none
 Args    : none


=cut

sub end_document{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 start_element

 Title   : start_element
 Usage   : $eventgenerator->start_element
 Function: Handles a start element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub start_element{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 end_element

 Title   : start_element
 Usage   : $eventgenerator->end_element
 Function: Handles an end element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub end_element{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}


=head2 in_element

 Title   : in_element
 Usage   : if( $eventgenerator->in_element($element) ) {}
 Function: Test if we are in a particular element
           This is different than 'within' because 'in' tests only
           if one has reached a specific element.
 Returns : boolean
 Args    : string element name 


=cut

sub in_element{
   my ($self,@args) = @_;
   $self->throw_not_implemented;

}

=head2 within_element

 Title   : within_element
 Usage   : if( $eventgenerator->within_element($element) ) {}
 Function: Test if we are within a particular element
           This is different than 'in' because within can be tested
           for a whole block.
 Returns : boolean
 Args    : string element name 


=cut

sub within_element{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 characters

 Title   : characters
 Usage   : $eventgenerator->characters($str)
 Function: Send a character events
 Returns : none
 Args    : string


=cut

sub characters{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

1;
