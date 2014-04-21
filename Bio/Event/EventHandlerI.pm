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

=head1 Developer Notes

EventHandlerI implementations are used in the BioPerl IO systems to
decouple the task of tokenizing the input stream into data elements
and their attributes, which is format-specific, and the task of
collecting those elements and attributes into whatever is the result
of a parser, which is specific to the kind of result to be produced,
such as BioPerl objects, a tabular or array data structure, etc.

You can think of EventHandlerI-compliant parsers as faking a SAX XML
parser, making their input (typically a non-XML document) behave as if
it were XML. The overhead to do this can be quite substantial, at the
gain of not having to duplicate the parsing code in order to change
the parsing result, and not having to duplicate the logic of
instantiating objects between parsers for different formats that all
give rise to the same types of objects. This is perhaps best
illustrated by the Bio::SearchIO system, where many different formats
exist for sequence similarity and pairwise sequence alignment exist
that essentially all result in Bio::Search objects.

The method names and their invocation semantics follow their XML SAX
equivalents, see http://www.saxproject.org/apidoc/, especially the
org.xml.sax.ContentHandler interface.

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
 Usage   : $resultObj = $parser->start_document();
 Function: Receive notification of the beginning of a document (the
           input file of a parser). The parser will invoke this method
           only once, before any other event callbacks.

           Usually, a handler will reset any internal state structures
           when this method is called.

 Returns : none
 Args    : none


=cut

sub start_document{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 end_document

 Title   : end_document
 Usage   : $parser->end_document();
 Function: Receive notification of the end of a document (normally the
           input file of a parser). The parser will invoke this method
           only once, and it will be the last method invoked during
           the parse of the document. The parser shall not invoke this
           method until it has either abandoned parsing (because of an
           unrecoverable error) or reached the end of input.

           Unlike the XML SAX signature of this method, this method is
           expected to return the object representing the result of
           parsing the document.

 Returns : The object representing the result of parsing the input
           stream between the calls to start_document() and this method.
 Args    : none


=cut

sub end_document{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 start_element

 Title   : start_element
 Usage   : $parser->start_element

 Function: Receive notification of the beginning of an element. The
           Parser will invoke this method at the beginning of every
           element in the input stream; there will be a corresponding
           end_element() event for every start_element() event (even when
           the element is empty). All of the element's content will be
           reported, in order, before the corresponding end_element()
           event.

 Returns : none
 Args : A hashref with at least 2 keys: 'Data' and 'Name'. The value
        for 'Name' is expected to be the type of element being
        encountered; the understood values will depend on the IO
        parser to which this interface is being applied. Likewise, the
        value for 'Data' will be specific to event handler
        implementions, and the specific data chunking needs of input
        formats to be handled efficiently.


=cut

sub start_element{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 end_element

 Title   : end_element
 Usage   : $parser->end_element

 Function: Receive notification of the end of an element. The parser
           will invoke this method at the end of every element in the
           input stream; there will be a corresponding start_element()
           event for every end_element() event (even when the element
           is empty).

 Returns : none

 Args    : hashref with at least 2 keys, 'Data' and 'Name'. The semantics
           are the same as for start_element().


=cut

sub end_element{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}


=head2 in_element

 Title   : in_element
 Usage   : if( $handler->in_element($element) ) {}

 Function: Test if we are in a particular element. 

           Normally, in_element() will test for particular attributes,
           or nested elements, within a containing
           element. Conversely, the containing element can be queries
           with within_element(). The names understood as argument
           should be the same as the ones understood for the 'Name'
           key in start_element() and end_element().

           Typically, handler implementations will call this method
           from within the characters() method to determine the
           context of the data that were passed to characters().

 Returns : boolean 

 Args    : A string, the name of the element (normally an attribute name or nested sub-element name). 

=cut

sub in_element{
   my ($self,@args) = @_;
   $self->throw_not_implemented;

}

=head2 within_element

 Title   : within_element
 Usage   : if( $handler->within_element($element) ) {}

 Function: Test if we are within a particular kind of element. 

           Normally, the element type names understood as argument
           values will be for containing elements or data
           chunks. Conversely, in_element() can be used to test
           whether an attribute or nested element is the ccurrent
           context.

           Typically, a handler will call this method from within the
           characters() method to determine the context for the data
           that were passed to characters().

 Returns : boolean
 Args    : string element name 


=cut

sub within_element{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

=head2 characters

 Title   : characters
 Usage   : $parser->characters($str)
 Function: Receive notification of character data. The parser will
           call this method to report values of attributes, or larger
           data chunks, depending on the IO subsystem and event
           handler implementation. Values may be whitespace-padded
           even if the whitespace is insignificant for the format.

           The context of the character data being passed can be
           determined by calling the in_element() and within_element()
           methods.

 Returns : none
 Args    : string, the character data


=cut

sub characters{
   my ($self,@args) = @_;
   $self->throw_not_implemented;
}

1;
