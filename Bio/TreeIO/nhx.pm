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
#   my $empty = chr(20);
 
   # replace empty labels with a tag
#   s/\(,/\($empty,/ig;
#   s/,,/,$empty,/ig;
#   s/,,/,/ig;
#   s/,\)/,$empty\)/ig;
#   s/\"/\'/ig;

   my $chars = '';
   $self->_eventHandler->start_document;
   my $lastevent = '';
   my @ch = split(//, $_);
   for (my $i = 0 ; $i < @ch ; $i++) {
       my $ch = $ch[$i];
       if( $ch eq ';' ) { 	   
	   return $self->_eventHandler->end_document;
       } elsif ($ch eq '[') {
	   if ( length $chars ) {
	       if ( $lastevent eq ':' ) {
		   $self->_eventHandler->start_element( { Name => 'branch_length' } );
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( { Name => 'branch_length' });
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
           $chars = '';
	   $self->_eventHandler->start_element( { Name => 'nhx_tag' });
	   $lastevent = $ch;
       } elsif( $ch eq '(' ) {
	   $chars = '';
	   $self->_eventHandler->start_element( {'Name' => 'tree'} );
	   $lastevent = $ch;
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
		   $self->debug("id with no branchlength is $chars\n");
		   $self->_eventHandler->start_element( { 'Name' => 'node' } );
		   $self->_eventHandler->start_element( { 'Name' => 'id' } );
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( { 'Name' => 'id' } );
	       }
	       
	   } elsif ( $lastevent ne ']' ) {
	       $self->_eventHandler->start_element( {'Name' => 'node'} )
	   }
	   $self->_eventHandler->end_element( {'Name' => 'node'} );
	   $self->_eventHandler->end_element( {'Name' => 'tree'} );
	   $chars = '';
	   $lastevent = $ch;
       } elsif ( $ch eq ',' ) {
	   if( length $chars ) {
	       if( $lastevent eq ':' ) {
		   $self->_eventHandler->start_element( { 'Name' => 'branch_length'});
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( {'Name' => 'branch_length'});
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
	   $lastevent = $ch;
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
	   $lastevent = $ch;
       } elsif ( $ch eq '=' ) {
	   if ($self->_eventHandler->within_element('nhx_tag')) {
	       $self->_eventHandler->start_element( { Name => 'tag_name' } );
	       $self->_eventHandler->characters($chars);
	       $self->_eventHandler->end_element( { Name => 'tag_name' } );
	       $chars = '';
	       $lastevent = $ch;
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
	       $lastevent = $ch;
	   } else {
	       $chars .= $ch;
	   }
       } else { 	   
	   $chars .= $ch;
       }
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
       if($data[-1] !~ /\)$/ ) {
	   $data[0] = "(".$data[0];
	   $data[-1] .= ")";
       }
       $self->_print(join(',', @data), ";\n");
   }
   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return;
}

sub _write_tree_Helper {
    my ($node) = @_;
    return () if (!defined $node);

    my @data;
    
    foreach my $n ( $node->each_Descendent() ) {
	push @data, _write_tree_Helper($n);
    }

    if( @data > 1 ) {
	$data[0] = "(" . $data[0];
	$data[-1] .= ")";
	$data[-1] .= ":". $node->branch_length if $node->branch_length;
	$data[-1] .= '[' . join(":", "&&NHX", map { "$_=" . $node->nhx_tag($_) } keys %{$node->nhx_tag || {}}) . ']' if $node->can('nhx_tag');
    } else {
	push @data, $node->to_string;
    }
    return @data;
}


1;
