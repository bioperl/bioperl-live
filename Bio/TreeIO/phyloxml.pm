# $Id: phyloxml.pm 11507 2007-06-23 01:37:45Z jason $
#
# BioPerl module for Bio::TreeIO::phyloxml
#
# Cared for by Mira Han <mirhan@indiana.edu>
#
# Copyright Mira Han
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::phyloxml - TreeIO implementation for parsing 
    PhyloXML format.

=head1 SYNOPSIS

  # do not use this module directly
  use Bio::TreeIO;
  my $treeio = Bio::TreeIO->new(-format => 'phyloxml', -file => 'tree.dnd');
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

This module handles parsing and writing of phyloXML format.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted viax the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Mira Han

Email mirhan@indiana.edu

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::phyloxml;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Tree::NodePhyloXML;
use Bio::Event::EventGeneratorI;
use XML::SAX;
use Bio::TreeIO::PhyloXMLEventBuilder;


use base qw(Bio::TreeIO);

sub _initialize {
  my($self, %args) = @_;
  $args{-nodetype} ||= 'Bio::Tree::NodePhyloXML';
  $self->SUPER::_initialize(%args);
  $self->debug("Creating obj phyloxml\n");
  $self->attach_EventHandler(Bio::TreeIO::PhyloXMLEventBuilder->new(-verbose => $self->verbose(), %args));
  $self->{'_parser'} = XML::SAX::ParserFactory->parser('Handler' => $self->{'_handler'});
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : Bio::Tree::TreeI
 Args    : none


=cut

sub next_tree 
{
  my ($self) = @_;
  local $/ = ";\n";
  return unless $_ = $self->_readline;

  $self->debug("entry is $_\n");
  $self->{'_parser'}->parse_string($_);


  my $chars = '';
  $self->_eventHandler->start_document;
  my $tree = $self->_eventHandler->end_document($chars);
  return $tree;
}

=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Write a tree out to data stream in phyloxml format
 Returns : none
 Args    : Bio::Tree::TreeI object

=cut

sub write_tree{
}

sub _write_tree_Helper {
}


1;
