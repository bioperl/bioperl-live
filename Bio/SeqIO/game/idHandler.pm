# $Id$
#
# BioPerl module for Bio::SeqIO::game::idHandler
#
# Cared for by Brad Marshall <bradmars@yahoo.com>
#         
# Copyright Brad Marshall
#
# You may distribute this module under the same terms as perl itself
# _history
# June 25, 2000     written by Brad Marshall
#
# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::game::idHandler - GAME helper via PerlSAX helper.

=head1 SYNOPSIS

GAME helper for parsing new ID objects from GAME XML. Do not use directly

=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and 
other Bioperl modules. Send your comments and suggestions preferably 
to one of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org        - Bioperl list
  bioxml-dev@bioxml.org        - Technical discussion - Moderate volume
  bioxml-announce@bioxml.org   - General Announcements - Pretty dead
  http://www.bioxml.org/MailingLists/         - About the mailing lists

=head1 AUTHOR - Brad Marshall

Email: bradmars@yahoo.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# This template file is in the Public Domain.
# You may do anything you want with this file.
#

package Bio::SeqIO::game::idHandler;
use Bio::Root::RootI;

use vars qw{ $AUTOLOAD };
use strict;

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    # initialize ids
    $self->{'ids'} = [];

    return $self;
}

=head2 start_document

 Title   : start_document
 Usage   : $obj->start_document
 Function: PerlSAX method called when a new document is initialized
 Returns : nothing
 Args    : document name

=cut

# Basic PerlSAX
sub start_document            {
    my ($self, $document) = @_;
}

=head2 end_document

 Title   : end_document
 Usage   : $obj->end_document
 Function: PerlSAX method called when a document is finished for cleaning up
 Returns : list of ids seen
 Args    : document name

=cut

sub end_document              {
    my ($self, $document) = @_;
    return $self->{'ids'};
}

=head2 start_element

 Title   : start_element
 Usage   : $obj->start_element
 Function: PerlSAX method called when a new element is reached
 Returns : nothing
 Args    : element object

=cut

sub start_element             {
    my ($self, $element) = @_;

    if ($element->{'Name'} eq 'bx-seq:seq') {
	if ($element->{'Attributes'}->{'bx-seq:id'}) {
	    push @{$self->{'ids'}}, $element->{'Attributes'}->{'bx-seq:id'};
	} else {
	    if ($self->can('warn')) {
		$self->warn('WARNING: Attribute bx-seq:id is required on bx-seq:seq. Sequence will not be parsed.');
	    } else {
		warn('WARNING: Attribute bx-seq:id is required on bx-seq:seq. Sequence will not be parsed.');
	    }
	}
    }
    return 0;
}

=head2 end_element

 Title   : end_element
 Usage   : $obj->end_element
 Function: PerlSAX method called when an element is finished
 Returns : nothing
 Args    : element object

=cut

sub end_element               {
    my ($self, $element) = @_;

}

=head2 characters

 Title   : characters
 Usage   : $obj->end_element
 Function: PerlSAX method called when text between XML tags is reached
 Returns : nothing
 Args    : text

=cut

sub characters   {
    my ($self, $text) = @_;
}


=head2 AUTOLOAD

 Title   : AUTOLOAD
 Usage   : do not use directly
 Function: autoload handling of missing DESTROY method
 Returns : nothing
 Args    : text

=cut

# Others
sub AUTOLOAD {
    my $self = shift;

    my $method = $AUTOLOAD;
    $method =~ s/.*:://;
    return if $method eq 'DESTROY';

    print "UNRECOGNIZED $method\n";
}

1;

__END__
