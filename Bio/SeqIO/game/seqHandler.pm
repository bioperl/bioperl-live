# $Id$
#
# BioPerl module for Bio::SeqIO::game::seqHandler
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

Bio::SeqIO::game::seqHandler - GAME helper via PerlSAX helper.

=head1 SYNOPSIS

GAME helper for parsing new Sequence objects from GAME XML. Do not use directly

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

package Bio::SeqIO::game::seqHandler;
use vars qw{ $AUTOLOAD @ISA };

use XML::Handler::Subs;
use Bio::PrimarySeq;

@ISA = qw(XML::Handler::Subs);

sub new {
    my ($caller,$seq) = @_;
    my $class = ref($caller) || $caller;
    my $self = bless ( {
	string => '',
	seq  => $seq,
    }, $class);
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
    $self->{'in_current_seq'} = 'false';    
    $self->{'Names'} = [];
    $self->{'string'} = '';
}

=head2 end_document

 Title   : end_document
 Usage   : $obj->end_document
 Function: PerlSAX method called when a document is finished for cleaning up
 Returns : list of sequences seen
 Args    : document name

=cut

sub end_document              {
    my ($self, $document) = @_;
    delete $self->{'Names'};
    return new Bio::PrimarySeq( -seq => $self->{'residues'},
				-moltype => $self->{'moltype'},
				-id => $self->{'seq'},
				-accession => $self->{'accession'},
				-desc => $self->{'desc'},
				-length => $self->{'length'},
				);

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

    push @{$self->{'Names'}}, $element->{'Name'};
    $self->{'string'} = '';

    if ($element->{'Name'} eq 'bx-seq:seq') {
	if ($element->{'Attributes'}->{'bx-seq:id'} eq $self->{'seq'}) {
	    $self->{'in_current_seq'} = 'true';
	    $self->{'moltype'} = $element->{'Attributes'}->{'bx-seq:type'};
	    $self->{'length'} =  $element->{'Attributes'}->{'bx-seq:length'};
	} else {
	    #This is not the sequence we want to import, but that's ok
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

    if ($self->{'in_current_seq'} eq 'true') {      
	if ($self->in_element('bx-seq:residues')) {
	    while ($self->{'string'} =~ s/\s+//) {};
	    $self->{'residues'} = $self->{'string'};
	}


	if ($self->in_element('bx-seq:name')) {
	    $self->{'string'} =~ s/^\s+//g;
	    $self->{'string'} =~ s/\s+$//;
	    $self->{'string'} =~ s/\n//g;
	    $self->{'name'} = $self->{'string'};
	}


	if ($self->in_element('bx-link:id')  && $self->within_element('bx-link:dbxref')) {
	    $self->{'string'} =~ s/^\s+//g;
	    $self->{'string'} =~ s/\s+$//;
	    $self->{'string'} =~ s/\n//g;
	    $self->{'accession'} = $self->{'string'};
	}

	if ($self->in_element('bx-seq:description')) {
	    $self->{'desc'} = $self->{'string'};
	}

	if ($self->in_element('bx-seq:seq')) {
	    $self->{'in_current_seq'} = 'false';
	}
    }

    pop @{$self->{'Names'}};

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
    $self->{'string'} .= $text->{'Data'};
}

=head2 in_element

 Title   : in_element
 Usage   : $obj->in_element
 Function: PerlSAX method called to test if state is in a specific element
 Returns : boolean
 Args    : name of element

=cut

sub in_element {
    my ($self, $name) = @_;

    return ($self->{'Names'}[-1] eq $name);
}

=head2 within_element

 Title   : within_element
 Usage   : $obj->within_element
 Function: PerlSAX method called to list depth within specific element
 Returns : boolean
 Args    : name of element

=cut

sub within_element {
    my ($self, $name) = @_;

    my $count = 0;
    foreach my $el_name (@{$self->{'Names'}}) {
	$count ++ if ($el_name eq $name);
    }

    return $count;
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
