#
# BioPerl module for Bio::TreeIO::tabtree
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::tabtree - A simple output format which displays a tree as an ASCII drawing

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = Bio::TreeIO->new(-file => 'input', -format => 'newick');
  my $out = Bio::TreeIO->new(-file => '>output', -format => 'tabtree');

  while( my $tree = $in->next_tree ) {
      $out->write_tree($tree);
  }

=head1 DESCRIPTION

This is a made up format just for outputting trees as an ASCII drawing.

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::tabtree;
use strict;

# Object preamble - inherits from Bio::Root::Root



use base qw(Bio::TreeIO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO::tabtree->new();
 Function: Builds a new Bio::TreeIO::tabtree object 
 Returns : Bio::TreeIO::tabtree
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
   my $line = _write_tree_Helper($tree->get_root_node,"");
   $self->_print($line. "\n");   
   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return;
}

sub _write_tree_Helper {
    my ($node,$indent) = @_;
    return unless defined $node;

    my @d = $node->each_Descendent();
    my $line = "";
    my ($i,$lastchild) = (0,scalar @d - 1);
    for my $n ( @d ) {
	if( $n->is_Leaf ) {
	    $line .= sprintf("%s| \n%s\\-%s\n",
			     $indent,$indent,$n->id || '');
	} else { 
	    $line .= sprintf("$indent|  %s\n",( $n->id ? 
					       sprintf("(%s)",$n->id) : ''));
	}
	my $new_indent = $indent . (($i == $lastchild) ? "| " : "  ");
	if( $n != $node ) {
	    # avoid the unlikely case of cycles
	    $line .= _write_tree_Helper($n,$new_indent);	
	}
    }
    return $line;
}

=head2 next_tree

 Title   : next_tree
 Usage   : 
 Function: Sorry not possible with this format
 Returns : none
 Args    : none


=cut

sub next_tree{
    $_[0]->throw("Sorry the format 'tabtree' can only be used as an output format at this time");
}

1;
