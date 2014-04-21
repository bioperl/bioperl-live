#
# BioPerl module for Bio::TreeIO::nhx
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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
  my $treeio = Bio::TreeIO->new(-format => 'nhx', -file => 'tree.dnd');
  my $tree = $treeio->next_tree;

=head1 DESCRIPTION

This module handles parsing and writing of Newick/New Hampshire eXtended (NHX) format.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted viax the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Aaron Mackey

Email amackey-at-virginia.edu

=head1 CONTRIBUTORS

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::nhx;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Tree::NodeNHX;
use Bio::Event::EventGeneratorI;
#use XML::Handler::Subs;

use base qw(Bio::TreeIO::newick);

sub _initialize {
  my($self, %args) = @_;
  $args{-nodetype} ||= 'Bio::Tree::NodeNHX';
  $self->SUPER::_initialize(%args);
}

sub _node_as_string {
  my $self = shift;
  my $node = shift;
  my $params = shift;
  
  my $label_stringbuffer = $self->SUPER::_node_as_string($node,$params);

  my @tags = $node->get_all_tags;
  if( scalar(@tags) > 0 ) {
    @tags = sort @tags;
    $label_stringbuffer .= '[' . 
      join(":", "&&NHX",
           map { "$_=" .join(',',$node->get_tag_values($_)) } 
           @tags ) . ']';
  }
  return $label_stringbuffer;
}

1;
