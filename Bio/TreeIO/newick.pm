# $Id$
#
# BioPerl module for Bio::TreeIO::newick
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::newick - TreeIO implementation for parsing 
  Newick/New Hampshire/PHYLIP format.

=head1 SYNOPSIS

  # do not use this module directly
  use Bio::TreeIO;
  my $treeio = new Bio::TreeIO(-format => 'newick', 
                               -file => 't/data/LOAD_Ccd1.dnd');
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

This module handles parsing and writing of Newick/PHYLIP/New Hampshire format.

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


package Bio::TreeIO::newick;
use vars qw(@ISA);
use strict;

use Bio::TreeIO;
use Bio::Event::EventGeneratorI;


@ISA = qw(Bio::TreeIO );


=head2 _initialize

 Title   : _initialize
 Args    : -print_count => boolean  default is false


=cut

sub _initialize { 
    my $self = shift;
    $self->SUPER::_initialize(@_);
    my ($print_count) = $self->_rearrange([qw(PRINT_COUNT)],
					  @_);
    $self->print_tree_count($print_count || 0);
}


=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : L<Bio::Tree::TreeI>
 Args    : none


=cut

sub next_tree{
   my ($self) = @_;
   local $/ = ";\n";
   return unless $_ = $self->_readline;
   s/[\r\n]//gs;
   my $despace = sub {my $dirty = shift; $dirty =~ s/\s+//gs; return $dirty};
   my $dequote = sub {my $dirty = shift; $dirty =~ s/^"?\s*(.+?)\s*"?$/$1/; return $dirty};
   s/([^"]*)(".+?")([^"]*)/$despace->($1) . $dequote->($2) . $despace->($3)/egsx;

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
   my ($prev_event,$lastevent,$id) = ('','','');
   foreach my $ch ( split(//,$_) ) {
       if( $ch eq ';' ) {
	   return $self->_eventHandler->end_document;
       } elsif( $ch eq '(' ) {
	   $chars = '';
	   $self->_eventHandler->start_element( {'Name' => 'tree'} );
       } elsif($ch eq ')' ) {
	   if( length($chars) ) {
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
		   $id = $chars;
	       }
	       my $leafstatus = 0;
	       if( $lastevent ne ')' ) {
		   $leafstatus = 1;
	       }

	       $self->_eventHandler->start_element({'Name' => 'leaf'});
	       $self->_eventHandler->characters($leafstatus);
	       $self->_eventHandler->end_element({'Name' => 'leaf'});
	       $id = '';
	   } else {
	       $self->_eventHandler->start_element( {'Name' => 'node'} );
	   }

	   $self->_eventHandler->end_element( {'Name' => 'node'} );
	   $self->_eventHandler->end_element( {'Name' => 'tree'} );
	   $chars = '';
       } elsif ( $ch eq ',' ) {
	   if( length($chars) ) {
	       if( $lastevent eq ':' ) {
		   $self->_eventHandler->start_element( { 'Name' => 'branch_length'});
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( {'Name' => 'branch_length'});
		   $lastevent = $prev_event;
		   $chars = '';		   
	       } else { 
		   $self->debug("id with no branchlength is $chars\n");
		   $self->_eventHandler->start_element( { 'Name' => 'node' } );
		   $self->_eventHandler->start_element( { 'Name' => 'id' } );
		   $self->_eventHandler->characters($chars);
		   $self->_eventHandler->end_element( { 'Name' => 'id' } );
		   $id = $chars;
	       }
	   } else {
	       $self->_eventHandler->start_element( { 'Name' => 'node' } );
	   }
	   my $leafstatus = 0;
	   if( $lastevent ne ')' ) {
	       $leafstatus = 1;
	   }
	   $self->_eventHandler->start_element({'Name' => 'leaf'});
	   $self->_eventHandler->characters($leafstatus);
	   $self->_eventHandler->end_element({'Name' => 'leaf'});
	   $self->_eventHandler->end_element( {'Name' => 'node'} );
	   $chars = '';
	   $id    = '';
       } elsif( $ch eq ':' ) {
	   $self->debug("id with a branchlength coming is $chars\n");
	   $self->_eventHandler->start_element( { 'Name' => 'node' } );
	   $self->_eventHandler->start_element( { 'Name' => 'id' } );	   
	   $self->_eventHandler->characters($chars);
	   $self->_eventHandler->end_element( { 'Name' => 'id' } );	   
	   $id = $chars;
	   $chars = '';
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
 Function: Write a tree out to data stream in newick/phylip format
 Returns : none
 Args    : L<Bio::Tree::TreeI> object

=cut

sub write_tree{
   my ($self,@trees) = @_;      
   if( $self->print_tree_count ){ 
       $self->_print(sprintf(" %d\n",scalar @trees));
   }
   foreach my $tree( @trees ) {
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
	# let's explicitly write out the bootstrap if we've got it
	my $b;
	if( defined ($b = $node->bootstrap) ) {
	    $data[-1] .= $b;
	} elsif( defined ($b = $node->id) ) {
	    $data[-1] .= $b;
	}
	$data[-1] .= ":". $node->branch_length if( defined $node->branch_length);
	
    } else {
	if( defined $node->id || defined $node->branch_length ) { 
	    push @data, sprintf("%s%s",
				defined $node->id ? $node->id : '', 
				defined $node->branch_length ? ":" .
				$node->branch_length : '');
	}
    }
    return @data;
}

=head2 print_tree_count

 Title   : print_tree_count
 Usage   : $obj->print_tree_count($newval)
 Function: Get/Set flag for printing out the tree count (paml,protml way)
 Returns : value of print_tree_count (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub print_tree_count{
    my $self = shift;
    return $self->{'_print_tree_count'} = shift if @_;
    return $self->{'_print_tree_count'} || 0;
}



1;
