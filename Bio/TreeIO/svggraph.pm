#
# BioPerl module for Bio::TreeIO::svg-graph
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Brian O'Connor
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::svggraph - A simple output format that converts a Tree object to an SVG output

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = new Bio::TreeIO(-file => 'input', -format => 'newick');
  my $out = new Bio::TreeIO(-file => '>output', -format => 'svggraph');

  while( my $tree = $in->next_tree ) {
      my $svg_xml = $out->write_tree($tree);
  }

=head1 DESCRIPTION

This outputs a tree as an SVG graphic using the SVG::Graph API

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Brian OConnor

Email brian.oconnor@excite.com

Describe contact details here

=head1 CONTRIBUTORS

Allen Day

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::svggraph;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::TreeIO;
use SVG::Graph;
use SVG::Graph::Data;
use SVG::Graph::Data::Tree;
use SVG::Graph::Data::Node;
use Bio::Tree::TreeI;
use Bio::Tree::Node;
use Tree::DAG_Node;


@ISA = qw(Bio::TreeIO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::TreeIO::svggraph();
 Function: Builds a new Bio::TreeIO::svggraph object 
 Returns : Bio::TreeIO::svggraph
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

}

=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Write a tree out to data stream in newick/phylip format
 Returns : none
 Args    : Bio::Tree::TreeI object

=cut

sub write_tree{
   my ($self,$tree) = @_;
   my $line = _write_tree_Helper($tree->get_root_node);
   $self->_print($line. "\n");
   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return;
}

sub _write_tree_Helper {
   my ($node) = @_;

   #this needs to be parameterized
   my $graph = SVG::Graph->new(width=>1600,height=>1000,margin=>30);

   my $group0 = $graph->add_frame;
   my $tree = SVG::Graph::Data::Tree->new;
   my $root = SVG::Graph::Data::Node->new;
   $root->name($node->id);
   _decorateRoot($root, $node->each_Descendent());
   $tree->root($root);
   $group0->add_data($tree);

   #this needs to be parameterized
   $group0->add_glyph('tree', stroke=>'black','stroke-width'=>2,'font-size'=>'10px');

   return($graph->draw);
}


=head2 decorateRoot

 Title   : _decorateRoot
 Usage   : internal methods
 Function:
 Example :
 Returns :
 Args    :


=cut

sub _decorateRoot{
  my $previousNode = shift;
  my @children = @_;
   foreach my $child (@children)
	 {
	   my $currNode = SVG::Graph::Data::Node->new;
	   $currNode->branch_label($child->id);
	   $currNode->branch_length($child->branch_length);
	   $previousNode->add_daughter($currNode);
	   _decorateRoot($currNode, $child->each_Descendent());
	 }
}

=head2 next_tree

 Title   : next_tree
 Usage   : 
 Function: Sorry not possible with this format
 Returns : none
 Args    : none


=cut

sub next_tree{
    $_[0]->throw("Sorry the format 'svggraph' can only be used as an output format");
}

1;
