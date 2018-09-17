#
# BioPerl module for Bio::SearchIO::EventHandlerI
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

Bio::SearchIO::EventHandlerI - An abstract Event Handler for Search Result parsing

=head1 SYNOPSIS

# do not use this object directly it is an interface
# See Bio::SearchIO::SearchResultEventBuilder for an implementation

    use Bio::SearchIO::SearchResultEventBuilder;
    my $handler = Bio::SearchIO::SearchResultEventBuilder->new();

=head1 DESCRIPTION

This interface describes the basic methods needed to handle Events
thrown from parsing a Search Result such as FASTA, BLAST, or HMMer.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::EventHandlerI;
use strict;
use Carp;


use base qw(Bio::Event::EventHandlerI);

=head2 start_result

 Title   : start_result
 Usage   : $handler->start_result($data)
 Function: Begins a result event cycle
 Returns : none 
 Args    : Type of Result

=cut

sub start_result {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 end_result

 Title   : end_result
 Usage   : $handler->end_result($data)
 Function: Ends a result event cycle
 Returns : Bio::Search::Result::ResultI object
 Args    : none


=cut

sub end_result{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 start_hsp

 Title   : start_hsp
 Usage   : $handler->start_hsp($data)
 Function: Start a HSP event cycle
 Returns : none
 Args    : type of element
           associated hashref

=cut

sub start_hsp{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 end_hsp

 Title   : end_hsp
 Usage   : $handler->end_hsp()
 Function: Ends a HSP event cycle
 Returns : Bio::Search::HSP::HSPI object
 Args    : type of event and associated hashref

=cut

sub end_hsp{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 start_hit

 Title   : start_hit
 Usage   : $handler->start_hit()
 Function: Starts a Hit event cycle
 Returns : none
 Args    : type of event and associated hashref


=cut

sub start_hit {
   my ($self,@args) = @_;
   $self->throw_not_implemented
}

=head2 end_hit

 Title   : end_hit
 Usage   : $handler->end_hit()
 Function: Ends a Hit event cycle
 Returns : Bio::Search::Hit::HitI object
 Args    : type of event and associated hashref


=cut

sub end_hit {
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 start_iteration

 Title   : start_iteration
 Usage   : $handler->start_iteration()
 Function: Starts an Iteration event cycle
 Returns : none
 Args    : type of event and associated hashref


=cut

sub start_iteration {
   my ($self,@args) = @_;
   $self->throw_not_implemented
}

=head2 end_iteration

 Title   : end_iteration
 Usage   : $handler->end_iteration()
 Function: Ends an Iterationevent cycle
 Returns : Bio::Search::Iteration::IterationI object
 Args    : type of event and associated hashref


=cut

sub end_iteration {
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 register_factory

 Title   : register_factory
 Usage   : $handler->register_factory('TYPE',$factory);
 Function: Register a specific factory for a object type class
 Returns : none
 Args    : string representing the class and
           Bio::Factory::ObjectFactoryI

See L<Bio::Factory::ObjectFactoryI> for more information

=cut

sub register_factory{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


=head2 factory

 Title   : factory
 Usage   : my $f = $handler->factory('TYPE');
 Function: Retrieves the associated factory for requested 'TYPE'
 Returns : a Bio::Factory::ObjectFactoryI
 Throws  : Bio::Root::BadParameter if none registered for the supplied type
 Args    : name of factory class to retrieve

See L<Bio::Factory::ObjectFactoryI> for more information

=cut

sub factory{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


=head2 Bio::Event::EventHandlerI methods

=cut

=head2 will_handle

 Title   : will_handle
 Usage   : if( $handler->will_handle($event_type) ) { ... }
 Function: Tests if this event builder knows how to process a specific event
 Returns : boolean
 Args    : event type name


=cut

=head2 SAX methods

See L<Bio::Event::EventHandlerI> for the additional SAX methods.

=cut


1;
