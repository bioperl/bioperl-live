# $Id$
#
# BioPerl module for Bio::SeqIO::game::featureHandler
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

Bio::SeqIO::game::featureHandler - GAME helper via PerlSAX helper.

=head1 SYNOPSIS

GAME helper for parsing new Feature objects from GAME XML. Do not use directly.

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

package Bio::SeqIO::game::featureHandler;

use Bio::SeqFeature::Generic;
use XML::Handler::Subs;

use vars qw{ $AUTOLOAD @ISA };
use strict;

@ISA = qw(XML::Handler::Subs);

sub new {
    my ($caller,$seq,$length,$type) = @_;
    my $class = ref($caller) || $caller;
    my $self = bless ({
	seq      => $seq,
	type     => $type,
	length   => $length,
	string   => '',
	feat     => {},
	feats    => [],
	comp_id  => 1,
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

    $self->{'Names'} = [];
    $self->{'Nodes'} = [];
    $self->{'feats'} = [];

}

=head2 end_document

 Title   : end_document
 Usage   : $obj->end_document
 Function: PerlSAX method called when a document is finished for cleaning up
 Returns : list of features seen
 Args    : document name

=cut

sub end_document              {
    my ($self, $document) = @_;

    delete $self->{'Names'};
    return $self->{'feats'};
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

    if ($self->in_element('bx-feature:seq_relationship')) {
	if (defined $element->{'Attributes'}->{'bx-feature:seq'} && 
	    defined $self->{'seq'} &&
	    $element->{'Attributes'}->{'bx-feature:seq'} eq $self->{'seq'}) {
	    $self->{'in_current_seq'} = 'true';
	} 
    }


    if ($self->in_element('bx-computation:computation')) {
	$self->{'feat'} = {};
	if (defined $element->{'Attributes'}->{'bx-computation:id'}) {
	    $self->{'feat'}->{'computation_id'} = $element->{'Attributes'}->{'bx-computation:id'};
	}  else {
	    $self->{'feat'}->{'computation_id'} = $self->{'comp_id'};
	    $self->{'comp_id'}++;
	}
    }

    if ($self->in_element('bx-feature:feature')) {
	if (defined $element->{'Attributes'}->{'bx-feature:id'}) {
	    $self->{'feat'}->{'id'} = $element->{'Attributes'}->{'bx-feature:id'};
	}
    }

    if ($self->in_element('bx-annotation:annotation')) {
	$self->{'feat'} = {};
	$self->{'feat'}->{'annotation_id'} = $element->{'Attributes'}->{'bx-annotation:id'};
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

    if ($self->in_element('bx-computation:program')) {
	$self->{'string'} =~ s/^\s+//g;
	$self->{'string'} =~ s/\s+$//;
	$self->{'string'} =~ s/\n//g;
	$self->{'feat'}->{'source_tag'} = $self->{'string'};
    }

    if ($self->in_element('bx-annotation:author')) {
	$self->{'string'} =~ s/^\s+//g;
	$self->{'string'} =~ s/\s+$//;
	$self->{'string'} =~ s/\n//g;
	$self->{'feat'}->{'source_tag'} = "Annotated by $self->{'string'}.";
    }

    if ($self->in_element('bx-feature:type')) {
	$self->{'string'} =~ s/^\s+//g;
	$self->{'string'} =~ s/\s+$//;
	$self->{'string'} =~ s/\n//g;
	$self->{'feat'}->{'primary_tag'} = $self->{'string'};
    }

    if ($self->in_element('bx-feature:start')) {
	$self->{'string'} =~ s/^\s+//g;
	$self->{'string'} =~ s/\s+$//;
	$self->{'string'} =~ s/\n//g;
	$self->{'feat'}->{'start'} = $self->{'string'};
    }

    if ($self->in_element('bx-feature:end')) {
	$self->{'string'} =~ s/^\s+//g;
	$self->{'string'} =~ s/\s+$//;
	$self->{'string'} =~ s/\n//g;
	$self->{'feat'}->{'end'} = $self->{'string'};
    }

    if ($self->in_element('bx-computation:score')) {
	$self->{'string'} =~ s/^\s+//g;
	$self->{'string'} =~ s/\s+$//;
	$self->{'string'} =~ s/\n//g;
	$self->{'feat'}->{'score'} = $self->{'string'};
    }

    if ($self->in_element('bx-feature:seq_relationship')) {
	if ($self->{'feat'}->{'start'} > $self->{'feat'}->{'end'}) {
	    my $new_start = $self->{'feat'}->{'end'};
	    $self->{'feat'}->{'end'} = $self->{'feat'}->{'start'};
	    $self->{'feat'}->{'start'} = $new_start;
	    $self->{'feat'}->{'strand'} = -1;
	} else {
	    $self->{'feat'}->{'strand'} = 1;
	}
	my $new_feat = new Bio::SeqFeature::Generic(
						    -start=>$self->{'feat'}->{'start'},
						    -end=>$self->{'feat'}->{'end'},
						    -source=>$self->{'feat'}->{'source_tag'},
						    -primary=>$self->{'feat'}->{'primary_tag'},
						    -score=>$self->{'feat'}->{'score'}
						    );
	if (defined $self->{'feat'}->{'computation_id'}) {
	    $new_feat->add_tag_value('computation_id', $self->{'feat'}->{'computation_id'} );
	} elsif (defined $self->{'feat'}->{'annotation_id'}) {
	    $new_feat->add_tag_value('annotation_id', $self->{'feat'}->{'annotation_id'} );
	}
	if (defined $self->{'feat'}->{'id'}) {
	    $new_feat->add_tag_value('id', $self->{'feat'}->{'id'} );
	}

	push @{$self->{'feats'}}, $new_feat;
	$self->{'feat'} = { 
	    seqid => $self->{'feat'}->{'curr_seqid'},
	    primary_tag => $self->{'feat'}->{'primary_tag'},
	    source_tag => $self->{'feat'}->{'source_tag'},
	    computation_id => $self->{'feat'}->{'computation_id'},
	    annotation_id => $self->{'feat'}->{'annotation_id'}
	}
    }


    pop @{$self->{'Names'}};
    pop @{$self->{'Nodes'}};

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
    
    return (defined $self->{'Names'}[-1] && 
	    $self->{'Names'}[-1] eq $name);
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
