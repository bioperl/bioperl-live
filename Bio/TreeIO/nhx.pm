# $Id$
#
# BioPerl module for Bio::TreeIO::nhx
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::nhx - TreeIO implementation for parsing 
    Newick/New Hampshire eXtendend (NHX) format.

=head1 SYNOPSIS

  # do not use this module directly
  use Bio::TreeIO;
  my $treeio = new Bio::TreeIO(-format => 'nhx', -file => 'tree.dnd');
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

This module handles parsing and writing of Newick/New Hampshire eXtended (NHX) format.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::nhx;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::TreeIO;
use Bio::Tree::NodeNHX;
use Bio::Event::EventGeneratorI;
#use XML::Handler::Subs;


@ISA = qw(Bio::TreeIO );

sub _initialize {
  my($self, %args) = @_;
  $args{-nodetype} ||= 'Bio::Tree::NodeNHX';
  $self->SUPER::_initialize(%args);
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : Bio::Tree::TreeI
 Args    : none


=cut

sub next_tree{
    my ($self) = @_;
    local $/ = ";\n";
    return unless $_ = $self->_readline;
    s/\s+//g;
    $self->debug("entry is $_\n");
    my $chars = '';
    $self->_eventHandler->start_document;
    my ($prev_event,$lastevent,$last_leaf_event) = ('','','');
    my @ch = split(//, $_);
    foreach my $ch  (@ch) {
	if( $ch eq ';' ) { 	   
	    $self->_eventHandler->in_element('node') && 
		$self->_eventHandler->end_element( {'Name' => 'node'});
	    return $self->_eventHandler->end_document;
	} elsif ($ch eq '[') {
	    if ( length $chars ) {
		if ( $lastevent eq ':' ) {
		    $self->_eventHandler->start_element( { Name => 'branch_length' } );
		    $self->_eventHandler->characters($chars);
		    $self->_eventHandler->end_element( { Name => 'branch_length' });
		    $lastevent = $prev_event;
		} else {
		    $self->debug("id with no branchlength is $chars\n");
		    $self->_eventHandler->start_element( { 'Name' => 'node' } );
		    $self->_eventHandler->start_element( { 'Name' => 'id' } );
		    $self->_eventHandler->characters($chars);
		    $self->_eventHandler->end_element( { 'Name' => 'id' } );		   
		}
	    } else {
		$self->_eventHandler->start_element( { Name => 'node' } );
	    }
	    my $leafstatus = ( $last_leaf_event ne ')' ) ? 1 : 0;
	    $self->_eventHandler->start_element({'Name' => 'leaf'});
	    $self->_eventHandler->characters($leafstatus);
	    $self->_eventHandler->end_element({'Name' => 'leaf'});	   
	    $chars = '';
	    
	    $self->_eventHandler->start_element( { Name => 'nhx_tag' });
	} elsif( $ch eq '(' ) {
	    $chars = '';
	    $self->_eventHandler->start_element( {'Name' => 'tree'} );
	} elsif($ch eq ')' ) {
	    if( length $chars ) {
		if( $lastevent eq ':') {
		    unless ($self->_eventHandler->within_element('nhx_tag')) {
			$self->_eventHandler->start_element( { 'Name' => 'branch_length'});
			$self->_eventHandler->characters($chars);
			$self->_eventHandler->end_element( {'Name' => 'branch_length'});
		    } else {
			$self->throw("malformed input; end of node ) before ] found");
		    }
		} else { 
		    $self->debug("id with no branchlength is '$chars'\n");
		    $self->_eventHandler->start_element( { 'Name' => 'node' } );
		    $self->_eventHandler->start_element( { 'Name' => 'id' } );
		    $self->_eventHandler->characters($chars);
		    $self->_eventHandler->end_element( { 'Name' => 'id' } );
		}
	    } elsif ( $lastevent ne ']' ) {
		$self->_eventHandler->start_element( {'Name' => 'node'} );
	    }
	    # problem here is that we need to detect if we coming up on
	    # the end of a leaf node or a labeled internal node
	    # each can have [] and each can have :, but only leaves are 
	    # NOT proceeded by a ')'
	    # the [] events throw us off
	    my $leafstatus = ( $last_leaf_event ne ')' ) ? 1 : 0;
	    $self->_eventHandler->start_element({'Name' => 'leaf'});
	    $self->_eventHandler->characters($leafstatus);
	    $self->_eventHandler->end_element({'Name' => 'leaf'});	   
	    
	    $self->_eventHandler->end_element( {'Name' => 'node'} );
	    $self->_eventHandler->end_element( {'Name' => 'tree'} );
	    $chars = '';
	    $last_leaf_event = $ch;

	} elsif ( $ch eq ',' ) {
	    if( length $chars ) {
		if( $lastevent eq ':' ) {
		    $self->_eventHandler->start_element( { 'Name' => 'branch_length'});
		    $self->_eventHandler->characters($chars);
		    $self->_eventHandler->end_element( {'Name' => 'branch_length'});
		    $lastevent = $prev_event;
		} else { 
		    $self->debug("id with no branchlength is $chars\n");
		    $self->_eventHandler->start_element( { 'Name' => 'node' } );
		    $self->_eventHandler->start_element( { 'Name' => 'id' } );
		    $self->_eventHandler->characters($chars);
		    $self->_eventHandler->end_element( { 'Name' => 'id' } );
		}   
	    } elsif ( $lastevent ne ']' ) {
		$self->_eventHandler->start_element( { 'Name' => 'node' } );
	    }
	    $self->_eventHandler->end_element( {'Name' => 'node'} );
	    $chars = '';
	    $last_leaf_event = $ch;
	} elsif( $ch eq ':' ) {
	    if ($self->_eventHandler->within_element('nhx_tag')) {
		if ($lastevent eq '=') {
		    $self->_eventHandler->start_element( { Name => 'tag_value' } );
		    $self->_eventHandler->characters($chars);
		    $self->_eventHandler->end_element( { Name => 'tag_value' } );
		    $chars = '';
		} else {
		    if ($chars eq '&&NHX') {
			$chars = ''; # get rid of &&NHX:
		    } else {
			$self->throw("Unrecognized, non \&\&NHX string: >>$chars<<");
		    }
		}
	    } elsif ($lastevent ne ']') {
		$self->debug("id with a branchlength coming is $chars\n");
		$self->_eventHandler->start_element( { 'Name' => 'node' } );
		$self->_eventHandler->start_element( { 'Name' => 'id' } );
		$self->_eventHandler->characters($chars);
		$self->_eventHandler->end_element( { 'Name' => 'id' } );
		$chars = '';
	    }
	} elsif ( $ch eq '=' ) {
	    if ($self->_eventHandler->within_element('nhx_tag')) {
		$self->_eventHandler->start_element( { Name => 'tag_name' } );
		$self->_eventHandler->characters($chars);
		$self->_eventHandler->end_element( { Name => 'tag_name' } );
		$chars = '';
	    } else {
		$chars .= $ch;
	    }
	} elsif ( $ch eq ']' ) {
	    if ($self->_eventHandler->within_element('nhx_tag') && $lastevent eq '=') {
		$self->_eventHandler->start_element( { Name => 'tag_value' } );
		$self->_eventHandler->characters($chars);
		$self->_eventHandler->end_element( { Name => 'tag_value' } );
		$chars = '';
		$self->_eventHandler->end_element( { Name => 'nhx_tag' } );
	    } else {
		$chars .= $ch;
		next;
	    }
	} else { 	   
	    $chars .= $ch;
	    next;
	}
	$prev_event = $lastevent;
	$lastevent = $ch;
    }       
    return undef;
}

=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Write a tree out to data stream in nhx format
 Returns : none
 Args    : Bio::Tree::TreeI object

=cut

sub write_tree{
    my ($self,@trees) = @_;
   foreach my $tree ( @trees ) {
       my @data = _write_tree_Helper($tree->get_root_node);
       # per bug # 1471 do not include enclosing brackets.
       # this is sort of cheating but it should work
       # remove first and last paren if the set ends in a paren
       if($data[-1] =~ s/\)$// ) {
	   $data[0] =~ s/^\(//;
       }
       $self->_print(join(',', @data), ";\n");
   }
   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return;
}

sub _write_tree_Helper {
    my ($node) = @_;
    return () unless defined $node;
    # rebless
    $node = bless $node,'Bio::Tree::NodeNHX';
    my @data;
    
    foreach my $n ( $node->each_Descendent() ) {
	push @data, _write_tree_Helper($n);
    }
    
    if( @data > 1 ) {
	$data[0] = "(" . $data[0];
	$data[-1] .= ")";
	$data[-1] .= ":". $node->branch_length if $node->branch_length;	
	# this is to not print out an empty NHX for the root node which is 
	# a convience for how we get a handle to the whole tree
	if( $node->ancestor || $node->id || defined $node->branch_length ) {
	    $data[-1] .= '[' . 
		join(":", "&&NHX", 
		     map { "$_=" .join(',',$node->get_tag_values($_)) } 
		     $node->get_all_tags() ) . ']';
	}
    } else { 
	push @data, $node->to_string; # a leaf
    }
    return @data;
}


1;
