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

Bio::TreeIO::newick - TreeIO implementation for parsing newick/New Hampshire format.

=head1 SYNOPSIS

# do not use this module directly
use Bio::TreeIO;
my $treeio = new Bio::TreeIO(-format => 'newick', -file => 'tree.dnd');
my $tree = $treeio->next_tree;

=head1 DESCRIPTION

Describe the object here

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


package Bio::TreeIO::newick;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::TreeIO;
use Bio::Event::EventGeneratorI;
use XML::Handler::Subs;


@ISA = qw(Bio::TreeIO );

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : Bio::Tree::TreeI
 Args    :


=cut

sub next_tree{
   my ($self,@args) = @_;
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
   foreach my $ch ( split(//,$_) ) {
       if( $ch eq ';' ) { 	   
	   return $self->_eventHandler->end_document;
       } elsif( $ch eq '(' ) {
	   $chars = '';
       } elsif($ch eq ')' || $ch eq ',' ) {
	   $self->_eventHandler->characters($chars);
	   $self->_eventHandler->end_element( {'Name' => 'branch_length'});
	   $self->_eventHandler->end_element( {'Name' => 'node'} );
	   $chars = '';
       } elsif( $ch eq ':' ) {
	   $self->_eventHandler->start_element( { 'Name' => 'node' } );
	   $self->_eventHandler->start_element( { 'Name' => 'id' } );	   
	   $self->_eventHandler->characters($chars);
	   $self->_eventHandler->end_element( { 'Name' => 'id' } );
	   $self->_eventHandler->start_element( { 'Name' => 'branch_length'});
	   $chars = '';
       } else { 	   
	   $chars .= $ch;
       } 
   }
   return undef;
}



1;
