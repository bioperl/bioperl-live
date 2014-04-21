#
# BioPerl module for Bio::TreeIO::svg-graph
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
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
  my $in = Bio::TreeIO->new(-file => 'input', -format => 'newick');
  my $out = Bio::TreeIO->new(-file => '>output', -format => 'svggraph');

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Brian OConnor

Email brian.oconnor-at-excite.com

=head1 CONTRIBUTORS

Allen Day
Guillaume Rousse, Guillaume-dot-Rousse-at-inria-dot-fr

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::svggraph;
use strict;

# Object preamble - inherits from Bio::Root::Root

use SVG::Graph;
use SVG::Graph::Data;
use SVG::Graph::Data::Tree;
use SVG::Graph::Data::Node;
use Bio::Tree::TreeI;
use Bio::Tree::Node;
use Tree::DAG_Node;


use base qw(Bio::TreeIO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO::svggraph->new();
 Function: Builds a new Bio::TreeIO::svggraph object
 Returns : Bio::TreeIO::svggraph
 Args    :-width    => image width (default 1600)
          -height   => image height (default 1000)
          -margin   => margin (default 30)
          -stroke   => stroke color (default 'black')
          -stroke_width=> stroke width (default 2)
          -font_size=> font size (default '10px')
          -nomalize => undef or 'log' (default is undef)

=cut

sub _initialize {
    my $self = shift;
    my ($width,$height,$margin,$stroke,
	$stroke_width,$font_size,
	$normalize) = $self->_rearrange([qw
					 (WIDTH
					  HEIGHT
					  MARGIN
					  STROKE
					  STROKE_WIDTH
					  FONT_SIZE
					  NORMALIZE)],
					@_);
    $self->{_width}        = $width || 1600;
    $self->{_height}       = $height || 1000;
    $self->{_margin}       = defined $margin ? $margin : 30;
    $self->{_stroke}       = $stroke || 'black';
    $self->{_stroke_width} = $stroke_width || 2;
    $self->{_font_size}    = $font_size || '10px';
    $self->{_normalize}    = $normalize || '';
    $self->SUPER::_initialize(@_);
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
   my $line = $self->_write_tree_Helper($tree->get_root_node);
   $self->_print($line. "\n");
   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return;
}

sub _write_tree_Helper {
   my ($self,$node) = @_;

   my $graph = SVG::Graph->new
       ('width'   => $self->{'_width'},
	'height'  => $self->{'_height'},
	'margin'  => $self->{'_margin'});

   my $group0 = $graph->add_frame;
   my $tree = SVG::Graph::Data::Tree->new;
   my $root = SVG::Graph::Data::Node->new;
   $root->name($node->id);
   $self->_decorateRoot($root, $node->each_Descendent());
   $tree->root($root);
   $group0->add_data($tree);

   $group0->add_glyph('tree',
		      'stroke'      =>$self->{'_stroke'},
		      'stroke-width'=>$self->{'_stroke_width'},
		      'font-size'   =>$self->{'_font_size'});

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

sub _decorateRoot {
    my ($self,$previousNode,@children) = @_;
    for my $child (@children) {
        my $currNode = SVG::Graph::Data::Node->new;

        # if no ID is set, the branch label is intentionally set blank (bug in SVG::Graph)
        my $id = $child->id || '';
        $currNode->branch_label($id);
        my $length = $child->branch_length;
        if ($self->{_normalize} eq 'log') {
            $length = log($length + 1);
        }

        $currNode->branch_length($length);
        $previousNode->add_daughter($currNode);
        $self->_decorateRoot($currNode, $child->each_Descendent());
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
